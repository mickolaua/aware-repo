"""
# ----------------------------------------------------------------------------- 
# Project:     AWARE
# Name:        aware.consumer.py 
# Purpose:     Alert message receiving and handling
# 
# Author:      npank 
# 
# Created:     2022-08-29 
# Copyright:   (c) 2004-2022 AWARE Developers
# ----------------------------------------------------------------------------
"""

from __future__ import annotations

import asyncio
import pickle
from contextlib import suppress
from datetime import datetime
from io import BytesIO
from threading import Lock, Thread
from typing import Any, Mapping, Optional, Sequence
from uuid import uuid4

import aiomisc
import gcn_kafka
import numpy as np
import pytz
from aiomisc import Service
from astropy import units as u
from astropy.time import Time
from confluent_kafka import (
    Consumer,
    KafkaError,
    KafkaException,
    Message,
    TopicPartition,
)
from confluent_kafka.error import ConsumeError
from sqlalchemy import inspect

from .. import alert, sql
from ..alert import AlertParsers
from ..alert.crossmatch import (
    crossmatch_alerts,
    crossmatch_alerts_by_name,
    replace_with_matched,
)
from ..alert.target_info import TargetInfo
from ..alert.util import add_retracted, is_retracted, max_id
from ..config import CfgOption
from ..credentials import Credentials
from ..data import AlertMessage, DataPackage, products_dir
from ..logger import log
from ..planning.main import max_area_trigger
from ..planning.planner import ObservationPlanner
from ..site import Telescopes, default_sites
from ..topic import TOPICS, full_topic_name_to_short, get_topics

# TODO: move topic option to aware.topic ?
topics = CfgOption(
    "topics",
    [t.name for t in get_topics()],
    list,
    comment="The list of alert topics to listen to",
)
timeout = CfgOption(
    "timeout",
    -1,
    float,
    comment="Timeout in seconds to wait for alert message consuming (default -1 means infinite timeout)",
)
start_date = CfgOption(
    "start_date",
    datetime.now(tz=pytz.UTC).isoformat(),
    datetime.fromisoformat,
    comment="The start date from which to listen for alerts",
)
max_tasks = CfgOption(
    "max_tasks",
    10,
    int,
    comment="Maximum number of asynchronous alert processing tasks",
)


@aiomisc.threaded
def parse_alert(alert_msg: bytes, topic: str) -> TargetInfo | None:
    """Parse the alert"""
    parser = AlertParsers.get(topic, None)

    if parser is not None:
        target_info = parser.parse_alert(alert_msg)
        return target_info
    else:
        return None


@aiomisc.threaded
def dump_alert_to_db(alert_msg: str, info: TargetInfo) -> None:
    # Add the alert to the database
    loc = info.localization
    try:
        ra, dec = loc.center().ra.deg, loc.center().dec.deg
        r = loc.error_radius().to_value("deg")
    except AttributeError:
        ra = dec = r = None

    loc = loc if loc else None

    alert_tab = sql.models.Alert(
        alert_message=alert_msg,
        ra_center=ra,
        dec_center=dec,
        error_radius=r,
        trigger_date=info.trigger_date,
        event=info.event,
        origin=info.origin,
        importance=info.importance,
        localization=pickle.dumps(loc),
    )

    with sql.models.dbconnect() as session:
        session.add(alert_tab)


def datetime_to_offset(dtime: datetime) -> int:
    """Convert datetime object to integer offset
    (in msec since the given date).

    Parameters
    ----------
    dtime : datetime
        a datetime object

    Returns
    -------
    offs: int
        an integer offset
    """
    return int(dtime.timestamp() * 1000)


def prepare_consumer(
    conf: dict[str, Any],
    credits: Credentials,
    topics: Sequence[str] = topics.value,
    start_offset: int = datetime_to_offset(start_date.value),
) -> Consumer:
    log.info("Reading messages from %s UT", start_date.value)
    log.info("Offset: %i", start_offset)
    consumer = gcn_kafka.Consumer(
        conf, client_id=credits.id, client_secret=credits.secret
    )
    partitions = [TopicPartition(topic, 0, start_offset) for topic in topics]
    offsets = consumer.offsets_for_times(partitions)
    consumer.assign(offsets)
    consumer.subscribe(topics)

    return consumer


@aiomisc.threaded
def poll_message(consumer: Consumer, timeout: float):
    return consumer.poll(timeout=timeout)


class ConsumeLoop(Service):
    """
    GCN Kafka asynchronous consumption loop.

    Parameters
    ----------
    consumer: Consumer
        a Kafka client that will be connected to the GCN for message
        consumption
    data_queue: Queue[DataPackage]
        a data queue for communication with Telegram Bot thread via DataPackage objects
    messages_count: int
        a count of messages that will be retrieved from Kafka in one iteration
    loop: AbstractEventLoop
        an abstract event loop
    """

    def __init__(
        self,
        consumer: Consumer,
        data_queue: asyncio.Queue[DataPackage] = asyncio.Queue(),
        messages_count: int = 10,
        loop: asyncio.AbstractEventLoop | None = None,
        max_tasks: int = 10,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.consumer = consumer
        self.data_queue = data_queue
        self.messages_count = messages_count or 1
        self.loop: asyncio.AbstractEventLoop = loop or asyncio.get_event_loop()
        self._lock = asyncio.Lock()
        self._max_commit_messages: int = 10
        self._commited_messages: int = 0
        self._max_tasks = max_tasks
        self._semaphore = asyncio.Semaphore(max_tasks)

    def _shutdown(self):
        """Shutdown the consume loop"""
        try:
            self._running = False
            self._on_shutdown()
        except Exception as e:
            log.exception(
                "\nerror occured on shutdown, see traceback for details", exc_info=e
            )

    def _on_shutdown(self):
        pending = asyncio.all_tasks()
        for task in pending:
            task.cancel()
            # Now we should await task to execute it's cancellation.
            # Cancelled task raises asyncio.CancelledError that we can
            # suppress:
            with suppress(asyncio.CancelledError):
                self.loop.run_until_complete(task)

    async def consume_messages(self):
        # ! Use recursion with sleep here, since we don't want block execution
        async with self._semaphore:
            while True:
                message = await self._poll_message(timeout=timeout.value)
                await self._process_message(message)
                await asyncio.sleep(0.0)

    async def start(self):
        """Start the consumer. This method is required by `Service`."""
        self._running = True
        await self.consume_messages()

    async def _poll_message(self, timeout: float) -> Message | None:
        """Consume a single message from Kafka.

        Parameters
        ----------
        consumer : str
            a consumer that will be consuming messages
        timeout : float
            a timeout for a message
        """
        try:
            message = await poll_message(self.consumer, timeout)

            log.debug(
                "ConsumeLoop._consume_message(): received message on %s",
                datetime.now().isoformat(),
            )
        except Exception as e:
            log.exception(
                "\nunable to consume the message, see traceback for details", exc_info=e
            )
            message = None

        return message

    def _validate_message(self, message: Message) -> bool:
        if not message or not message.value():
            log.exception("\nmessage is empty")
            return False

        if message.error() and message.error().code() == KafkaError._PARTITION_EOF:
            log.exception("\nreached the end of the partition")
            return False

        return True

    async def _process_message(self, message: Message):  # timeout_sec: float):
        # message = await self._poll_message(self.consumer, timeout_sec)

        if not self._validate_message(message):
            log.exception(
                "Message has not passed validation, processing will not occur"
            )
            return

        que = self.data_queue

        text = message.value().decode("utf-8")
        log.info("\n" + "#" * 8 + "\n" + "New message\n" + str(datetime.now()))
        log.debug("Message topic: %s\n content: \n%s", message.topic(), message.value())

        # Retrieve target info
        log.debug("\nParsing the alert ...")
        info = await parse_alert(message.value(), message.topic())
        if info is None:
            log.info(
                "Alert was ignored. Possibly, it is only available in develop mode."
            )
            return
        else:
            # Display sender information in log with trigger date
            log.info(
                "Alert from %s (%s) triggered at %s",
                info.origin,
                info.event,
                info.trigger_date,
            )

        log.debug(info)

        # Query alert database to determine if alert should be retracted
        if not info.rejected:
            info.rejected = is_retracted(info.event, info.origin)

        # Do not send plots and observational messages if the event matches previous
        # ones, but does not refines localization area
        send_obs_data = not info.rejected
        refines_localization = True

        # Store rejected (or retracted) event
        if info.rejected:
            log.debug("The senter rejected the event.")
            m_id = max_id()
            if not m_id:
                m_id = 1
            else:
                m_id += 1
            alert.util.add_retracted(m_id, info.event, info.origin)

        @aiomisc.threaded
        def spatial_xmatch(info: TargetInfo) -> list[dict]:
            try:
                matched_alerts = crossmatch_alerts(info)
            except Exception as e:
                log.error("Failed to crossmatch alert in space: %s", e)
                matched_alerts = []
            return matched_alerts

        @aiomisc.threaded
        def name_xmatch(info: TargetInfo) -> list[dict]:
            try:
                matched_alerts = crossmatch_alerts_by_name(info.origin, info.event)
            except Exception as e:
                log.error("Failed to crossmatch alert by name: %s", e)
                matched_alerts = []
            return matched_alerts

        if info.localization:
            matched_alerts = await spatial_xmatch(info)
        else:
            send_obs_data = False
            matched_alerts = await name_xmatch(info)

        event = info.event
        origin = info.origin

        refines_localization = False
        if matched_alerts:
            # sorted_alerts = sorted(matched_alerts, key=lambda _a: _a["trigger_date"])
            earliest_alert = min(matched_alerts, key=lambda _a: _a["trigger_date"])
            origin = earliest_alert["origin"]
            event = earliest_alert["event"]
            log.debug(
                "Event %s %s matches %d already stored events (firstly mentioned as "
                "%s %s)",
                info.origin,
                info.event,
                len(matched_alerts),
                origin,
                event,
            )

            log.debug("All cross-matched events will be rejected.")
            for al in matched_alerts:
                if info.rejected:
                    alert.util.add_retracted(al["id"], al["event"], al["origin"])
                if is_retracted(al["event"], al["origin"]):
                    info.rejected = True

            if not info.rejected:
                if info.localization:
                    matched_error_radii = [a["error_radius"] for a in matched_alerts]
                    i_min = np.argmin(matched_error_radii)
                    matched_error_radius = matched_error_radii[i_min]
                    matched_localization = matched_alerts[i_min]["localization"]
                    error_radius = info.localization.error_radius().to_value(u.deg)

                    if error_radius < matched_error_radius:
                        log.debug(
                            "%s %s refines localization radius for %d already stored"
                            "events",
                            origin,
                            event,
                            len(matched_alerts),
                        )

                        refines_localization = True

                        # If matched error radius is already small enough, no need to
                        # run the planner again. Since, skymap is wholy inside FOV of a
                        # telescope anyway.
                        SMALL_RADIUS = 1.5 / 60  # 1.5 arcmin
                        refines_small = matched_error_radius < SMALL_RADIUS

                        # The same, when newer sky map is fully contained inside old one
                        # borders, and its area not significantly smaller than the old
                        # one. We probably already have plans for that skymap region.
                        matched_moc = pickle.loads(matched_localization).moc()
                        new_moc = info.localization.moc()
                        SKY_MAP_AREA_FRACTION = 0.75
                        new_inside_old = (
                            new_moc.union(matched_moc) == new_moc
                            and new_moc.sky_fraction
                            > SKY_MAP_AREA_FRACTION * matched_moc.sky_fraction
                        )

                        if refines_small or new_inside_old:
                            send_obs_data = False
                        else:
                            send_obs_data = True

                    else:
                        send_obs_data = False
                        refines_localization = False

        # Send the alert message before making plots and planning
        uuid = uuid4().hex

        # Save it to the local DB
        log.debug("\nDumping the alert into the database ...")
        await dump_alert_to_db(text, info)

        # Create alert message for Telegram subscribers.
        # The alert is not presented in the database, we consider it new.
        # When the alert that matches at-least one alert written in the database
        # we provide first mention of the event origin and name.
        # If the alert provides more neat localiztion area we notify user on that and
        # send them the updated observational data in upcoming messages.
        # async with self._lock:
        senter = info.origin
        event_head = "Event" if matched_alerts else "New event"
        trigger_ = info.trigger_date.replace(tzinfo=None).isoformat("T", "seconds")
        pack_type = full_topic_name_to_short(info.packet_type)
        body = (
            f"[{event_head}]\n"
            + f"From: {info.origin}\n"
            + f"Trigger ID: {info.event}\n"
            + f"Trigger Date: {trigger_ } UTC\n"
            + f"Packet: {pack_type}\n"
            + (f"First mention: {origin} {event}\n" if event_head == "Event" else "")
            + (f"{info.description}\n" if info.description else "")
            + (
                (
                    f"Localization area for {origin} {event} refined."
                    if event_head == "Event"
                    and refines_localization
                    and info.localization
                    else ""
                )
                if not info.rejected
                else "EVENT CANCELED!"
            )
        )

        first_mention = f"{origin} {event}" if event_head == "Event" else ""
        alert_message = AlertMessage(
            id=uuid,
            senter=senter,
            body=body,
            alert_type=info.packet_type,
            first_mention=first_mention,
            is_new=not matched_alerts,
            update=refines_localization,
            target_info=info,
        )
        await que.put(alert_message)

        def is_plannable() -> bool:
            loc = getattr(info, "localization", None)
            is_passed = not info.rejected
            has_skymap = loc is not None
            max_area = max_area_trigger.value
            curr_area = getattr(loc, "area", None)
            is_skymap_size_in_bounds = (
                curr_area is not None and curr_area().to_value("deg2") < max_area
            )
            plannable = (
                send_obs_data and is_passed and has_skymap and is_skymap_size_in_bounds
            )
            return plannable

        if is_plannable():
            now = Time(datetime.now().isoformat(), format="isot")

            telescopes = {name: Telescopes[name] for name in default_sites.value}

            wide_field_telescopes = [
                t for t in telescopes.values() if t.fov.is_widefield
            ]

            # We do not want to observe GLADE+ galaxies, when no distance constraints
            # are given.
            if info.origin.lower()[:3] in {"lvk", "lvc"}:
                narrow_field_telescopes = [
                    t for t in telescopes.values() if not t.fov.is_widefield
                ]
            else:
                narrow_field_telescopes = []
            wide_field_comment = "Observe these sky fields in the given order."
            narrow_field_comment = "Observe these galaxies in the given order."
            single_scan_comment = "Aim the telescope at the localization center."
            planner = ObservationPlanner(info, working_directory=products_dir.value)

            # Run planner in the separate thread, because it is a blocking function
            @aiomisc.threaded
            def plan_func(planner: ObservationPlanner) -> dict:
                plan_dicts = {}
                try:
                    planner.plan_observations(
                        epoch=now,
                        wide_field_telescopes=wide_field_telescopes,
                        narrow_field_telescopes=narrow_field_telescopes,
                        day=1,
                        narrow_field_comment=narrow_field_comment,
                        wide_field_comment=wide_field_comment,
                        single_scan_comment=single_scan_comment,
                    )
                    plan_result = planner.save_or_update_planning()
                    json_plan_reg_fn, plan_data = plan_result
                    plan_dicts = plan_data.get("plans", [])
                except Exception as e:
                    log.error("Error while planning observations: %s", e)

                return plan_dicts

            # Plan only for not-empty localization skymaps
            if info.localization is not None:
                plan_dicts = await plan_func(planner)
            else:
                plan_dicts = {}

            # Send the plans to the communication thread (e.g. Telegram)
            if plan_dicts:
                for plan_day in plan_dicts:
                    if plan_day["day"] == 1:
                        try:
                            data_package = DataPackage(
                                id=uuid,
                                alert_type=info.packet_type,
                                site=telescopes[plan_day["site_id"]],
                                plan_filename=plan_day["file"],
                                plot_fname=plan_day["plot_file"],
                                comment=plan_day["comment"],
                                target_info=info,
                            )
                        except (LookupError, TypeError, ArithmeticError):
                            # If registry contains site_id not presented in
                            # default_sites
                            pass
                        else:
                            await que.put(data_package)

        self._commited_messages += 1
        if self._commited_messages % self._max_commit_messages:
            self.consumer.commit(asynchronous=True)
