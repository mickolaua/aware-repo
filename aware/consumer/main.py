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
import os
import pickle
import posixpath
import sys
from contextlib import suppress
from datetime import datetime
from hashlib import md5
from io import BytesIO, StringIO
from pickle import dumps
from queue import Queue
from threading import Lock, Thread
from typing import Any, Mapping, Optional, Sequence
from uuid import uuid4

import gcn_kafka
import numpy as np
import pytz
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from confluent_kafka import (
    Consumer,
    KafkaError,
    KafkaException,
    Message,
    TopicPartition,
)
from confluent_kafka.error import ConsumeError

# from aware.sftp import upload_files
from pathvalidate import sanitize_filename, sanitize_filepath
from sqlalchemy import inspect

from aware.field import Field

from .. import alert, site, sql
from ..alert import AlertParsers
from ..alert.crossmatch import (
    crossmatch_alerts,
    crossmatch_alerts_by_name,
    replace_with_matched,
)
from ..alert.target_info import TargetInfo
from ..config import CfgOption
from ..credentials import Credentials
from ..data import TelegramAlertMessage  # , TelegramSFTPUrl
from ..data import TelegramDataPackage
from ..glade import GladeCatalog, GladeGalaxy
from ..json import JSON
from ..logger import log
from ..planning.planner import SkymapPlanner
from ..planning.program import create_observation_program
from ..planning.main import max_area_trigger
from ..site import Site, Telescopes
from ..topic import TOPICS, full_topic_name_to_short
from ..io import create_event_folder
from ..data import products_dir

topics = CfgOption("topics", TOPICS, list)
timeout = CfgOption("timeout", -1, int)
start_date = CfgOption(
    "start_date", datetime.now(tz=pytz.UTC).isoformat(), datetime.fromisoformat
)
default_site = CfgOption("default_site", "mondy_azt33ik", str)
ntasks = CfgOption("ntasks", 8, int)
max_localization_radius = CfgOption("max_localization_radius", 1, float)


running = True


async def shutdown(loop: asyncio.BaseEventLoop):
    global running
    running = False
    loop.stop()
    loop.close()
    sys.exit(0)


async def save_file(filename: str, content: str):
    try:
        with open(filename, "w+") as f:
            f.write(content)
    except Exception as e:
        log.exception("unable to save file %s", filename, exc_info=e)
    else:
        log.info("filed saved: %s", filename)


async def parse_alert(alert_msg: bytes, topic: str) -> TargetInfo | None:
    """Parse the alert"""
    parser = AlertParsers.get(topic, None)

    if parser is not None:
        target_info = parser.parse_alert(alert_msg)
        return target_info
    else:
        return None


async def dump_alert_to_db(alert_msg: str, info: TargetInfo) -> None:
    # Create SQLite session
    engine, session = sql.models.create_session()

    # Create alert table if it is not exist (e.g. fresh database)
    ins = inspect(engine)
    if not ins.has_table("alert"):
        sql.models.Alert.__table__.create(bind=engine, checkfirst=True)

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

    with session:
        session.add(alert_tab)
        session.commit()


async def get_glade_galaxies(info: TargetInfo) -> list[GladeGalaxy] | None:
    center = SkyCoord(info.ra_center * u.deg, info.dec_center * u.deg)
    radius: u.Unit = info.error_radius * u.deg
    galaxies = GladeCatalog.query_field(center, radius)

    return galaxies


async def consume_message(consumer: str, timeout: float) -> Message | None:
    """Consume a single message in an asynchronous manner.

    Parameters
    ----------
    consumer : str
        a consumer that will be consuming messages
    timeout : float
        a timeout for a message
    """
    # Aquire a single message
    try:
        message = consumer.poll(timeout=timeout)
    except Exception as e:
        log.exception(e)
        message = None

    return message


async def sort_galaxies_for_obs(
    galaxies: Sequence[FixedTarget], site: Site
) -> list[FixedTarget] | None:
    if not galaxies:
        return []

    start, stop = await site.nearest_observation_window(Time(datetime.now()))
    obs_targets = await site.observable_targets(galaxies, start, stop)

    log.debug("Start: %s", start.datetime)
    log.debug("End: %s", stop.datetime)

    if obs_targets:
        ordered_targets = await site.observation_order(
            obs_targets, start_time=start, end_time=stop
        )
    else:
        ordered_targets = None

    return ordered_targets


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
    topics: Sequence[str] = topics.get_value(),
    start_offset: int = datetime_to_offset(start_date.get_value()),
) -> Consumer:
    log.info("Reading messages from %s UT", start_date.get_value())
    log.info("Offset: %i", start_offset)
    consumer = gcn_kafka.Consumer(
        conf, client_id=credits.id, client_secret=credits.secret
    )
    partitions = [TopicPartition(topic, 0, start_offset) for topic in topics]
    offsets = consumer.offsets_for_times(partitions)
    consumer.assign(offsets)
    consumer.subscribe(topics)

    return consumer


class ConsumeLoop:
    """
    GCN Kafka asynchronous consumption loop.

    Parameters
    ----------
    consumer: Consumer
        a Kafka client that will be connected to the GCN for message
        consumption
    data_queue: Queue[DataPackage]
        a data queue for communication with Telegram Bot thread via
        DataPackage objects
    ntasks: int
        a number of tasks that can be performed in parallel (not dynamically
        changeable)
    """

    def __init__(
        self,
        consumer: Consumer,
        data_queue: asyncio.Queue[TelegramDataPackage],
        ntasks: int = ntasks.get_value(),
        max_commit_messages: int = 10,
    ) -> None:
        self._consumer = consumer
        self._data_queue = data_queue
        self._running = False
        self._loop = asyncio.get_event_loop()
        self._semaphore = asyncio.Semaphore(ntasks)
        self._lock = asyncio.Lock()
        self._max_commit_messages = max_commit_messages
        self._commited_messages = 0

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
                self._loop.run_until_complete(task)

    async def _do_run(self):
        async with self._semaphore:
            asyncio.ensure_future(self._do_run(), loop=self._loop)
            await self._consume_messages()

    def run(self):
        try:
            asyncio.ensure_future(self._do_run(), loop=self._loop)
            self._loop.run_forever()
        except Exception as e:
            log.exception(
                "\ncan not run the consumption loop, see traceback for " + "details",
                exc_info=e,
            )
            self._shutdown()

    async def _consume_messages(self, timeout_sec: float = timeout.get_value()):
        self._running = True
        await self._process_message(timeout_sec)

    async def _consume_message(self, consumer: Consumer) -> Message | None:
        """Consume a single message from Kafka.

        Parameters
        ----------
        consumer : str
            a consumer that will be consuming messages
        timeout : float
            a timeout for a message
        """
        try:
            message = consumer.poll(timeout=timeout.get_value())
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

    async def _validate_message(self, message: Message) -> bool:
        if not message or not message.value():
            log.exception("\nmessage is empty")
            return False

        if message.error() and message.error().code() == KafkaError._PARTITION_EOF:
            log.exception("\nreached the end of the partition")
            return False

        return True

    async def _process_message(self, timeout_sec: float):
        message = await self._consume_message(self._consumer)

        if not await self._validate_message(message):
            log.exception(
                "Message has not passed validation, processing will not occur"
            )
            return

        que = self._data_queue

        text = message.value().decode("utf-8")
        log.info("\n" + "#" * 8 + "\n" + "New message\n" + str(datetime.now()))
        log.debug("Message topic: %s\n content: \n%s", message.topic(), message.value())

        # Retrieve target info
        log.debug("\nParsing the alert ...")
        info = await parse_alert(message.value(), message.topic())
        if info is None:
            log.debug(
                "Alert was ignored. Possibly, it is only available in develop mode."
            )
            return
        log.debug(info)

        # Do not send plots and observational messages if the event matches previous
        # ones, but does not refines localization area
        send_obs_data = not info.rejected
        refines_localization = True

        # Store rejected (or retracted) event
        if info.rejected:
            log.debug("The senter rejected the event.")
            m_id = alert.util.max_id()
            if not m_id:
                m_id = 1
            else:
                m_id += 1
            alert.util.add_retracted(m_id, info.event, info.origin)

        # Before sending to subscribers, check if this event matches older ones
        if info.localization is not None:
            matched_alerts = await crossmatch_alerts(info)
        else:
            send_obs_data = False
            matched_alerts = await crossmatch_alerts_by_name(info.origin, info.event)

        event = info.event
        origin = info.origin

        if matched_alerts:
            sorted_alerts = sorted(matched_alerts, key=lambda _a: _a.trigger_date)
            origin = sorted_alerts[0].origin
            event = sorted_alerts[0].event
            log.debug(
                "Event %s %s matches %d already stored events (firstly mentioned as "
                "%s %s)",
                info.origin,
                info.event,
                len(matched_alerts),
                origin,
                event,
            )
            if info.rejected:
                log.debug("All cross-matched events will be rejected.")
                for al in matched_alerts:
                    alert.util.add_retracted(al.id, al.event, al.origin)
            else:
                if info.localization:
                    matched_error_radii = [a.error_radius for a in matched_alerts]
                    i_min = np.argmin(matched_error_radii)
                    matched_error_radius = matched_error_radii[i_min]
                    error_radius = info.localization.error_radius().to_value(u.deg)

                    if error_radius < matched_error_radius:
                        log.debug(
                            "%s %s refines localization radius for %d already stored"
                            "events",
                            origin,
                            event,
                            len(matched_alerts),
                        )
                        send_obs_data = True
                        refines_localization = True
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
        async with self._lock:
            senter = info.origin
            event_head = "Event" if matched_alerts else "New event"
            trigger_ = info.trigger_date.replace(tzinfo=None).isoformat("T", "seconds")
            pack_type = full_topic_name_to_short(info.packet_type)
            body = (
                f"[{event_head}]\n"
                + f"From: {origin}\n"
                + f"Trigger ID: {info.event}\n"
                + f"Trigger Date: {trigger_ } UTC\n"
                + f"Packet: {pack_type}\n"
                + (
                    f"First mention: {origin} {event}\n"
                    if event_head == "Event"
                    else ""
                )
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
                    else "EVENT IS RETRACKED!"
                )
            )
            alert_message = TelegramAlertMessage(uuid, senter, body, info.packet_type)
            await que.put(alert_message)

        # We can send observational data
        if (
            send_obs_data
            and not info.rejected
            and (info.localization.area().value < max_area_trigger.value)
        ):
            sites = list(
                s
                for s in site.Telescopes.values()
                if s.name in site.default_sites.value
            )
            now = Time(datetime.now().isoformat(), format="isot")
            loc = info.localization
            local_files = []
            dist_files = []

            # TODO:
            # For those scopes that has large enough FOV to cover the entire sky map
            # send airmass plot, otherwise plan multiple scope observations

            # Temporaly plan multiple scopes on LVK sky maps only
            if origin in {"LVC", "LVK"}:
                planner = SkymapPlanner(
                    loc.data,
                    event_name=f"{origin}_{event}",
                    working_directory=create_event_folder(
                        origin, event, products_dir.value
                    ),
                )
                wide_field_telescopes = list(
                    s.name for s in sites if s.fov.is_widefield
                )
                narrow_field_telescopes = list(
                    s.name for s in sites if not s.fov.is_widefield
                )
                planner.plan_observations(
                    wide_field_telescopes, narrow_field_telescopes
                )
                planner.save_plan_fits(wide_field_telescopes, narrow_field_telescopes)
                saved_blocks = planner.save_blocks()

                all_blocks = planner.field_blocks
                all_blocks.update(planner.target_blocks)
                for name, fields in all_blocks.items():
                    if fields:
                        ax = loc.plot(
                            site.Telescopes[name],
                            start=None,
                            stop=None,
                            targets=fields,
                        )
                        outdir = os.path.join(planner.working_directory, name)
                        plot_name = (
                            f"plan_map_{origin}_{event}_{name}_day{planner.day}.png"
                        )
                        plot_name_safe = sanitize_filename(
                            plot_name, replacement_text="_"
                        )
                        plot_path = os.path.join(outdir, plot_name_safe)
                        fig = ax.get_figure()
                        fig.savefig(plot_path)

                        async with self._lock:
                            data_package = TelegramDataPackage(
                                uuid,
                                info,
                                info.packet_type,
                                site.Telescopes[name],
                                saved_blocks[name],
                                plot_path,
                            )
                            await que.put(data_package)

            else:
                for i, s in enumerate(sites):
                    start, stop = s.nearest_observation_window(now)

                    # Can not be observed due to its location
                    if start is None or stop is None:
                        continue

                    # If nearest observational window is located after 1 day since
                    # today, it is considered that targets are not observable
                    if start.jd - now.jd >= 1:
                        sorted_targets = []
                        comment = (
                            "Targets do not cross the horizon at "
                            f"{s.full_name} within 24 hours."
                        )
                    else:
                        sorted_targets, comment = await loc.observe(s, start, stop)

                    # Dump sorted targets
                    if sorted_targets:
                        obs_program = create_observation_program(s, sorted_targets)
                        ax = loc.plot(s, start, stop, sorted_targets)
                    else:
                        obs_program = ""

                    # If observational program is not empty
                    if obs_program:
                        outdir = os.path.join(
                            products_dir.get_value(),
                            f"{origin}_{event}",
                            f"{s.name}",
                        )

                        # Each event will be located in a separate subfolder
                        if not os.path.exists(outdir):
                            os.makedirs(outdir)

                        event_name = f"targets_{origin}_{event}_{s.name}"
                        fname = os.path.join(
                            outdir,
                            f"{event_name}.list",
                        )
                        fname_safe = sanitize_filepath(fname, replacement_text="_")
                        await save_file(fname_safe, obs_program)

                        plot_fname = os.path.join(outdir, f"{event_name}.png")
                        plot_fname_safe = sanitize_filepath(
                            plot_fname, replacement_text="_"
                        )
                        fig = ax.get_figure()
                        fig.savefig(plot_fname)
                    else:
                        # Do not send empty observational program
                        fname = ""
                        plot_fname = ""

                    async with self._lock:
                        data_package = TelegramDataPackage(
                            uuid, info, info.packet_type, s, fname, plot_fname
                        )
                        await que.put(data_package)

        #     if obs_program:
        #         upload_files(fname, f"/uploads/{info.event}/{s.name}/targets_{uuid}_{s.name}.{s.default_target_list_fmt}")
        #         upload_files(plot_fname, f"/uploads/{info.event}/{s.name}/targets_{uuid}_{s.name}.png")

        # async with self._lock:
        #     url = TelegramSFTPUrl(uuid, info.event, "", f"/uploads/{uuid}")
        #     await que.put(url)

        # async with self._lock:
        #     data_package = TelegramDataPackage(
        #         uuid, info, s.full_name, fname, plot_fname
        #     )
        #     await que.put(data_package)

        self._commited_messages += 1
        if self._commited_messages % self._max_commit_messages:
            self._consumer.commit(asynchronous=True)

    def to_thread(self) -> Thread:
        return Thread(name="consumer", target=self.run, daemon=True)
