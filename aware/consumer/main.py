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
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from confluent_kafka import (Consumer, KafkaError, KafkaException, Message,
                             TopicPartition)
from confluent_kafka.error import ConsumeError
from geojson import GeoJSON, dumps
from geojson.geometry import MultiPoint, Point
from matplotlib.axes import Axes

from aware.credentials import Credentials
from aware.data import TelegramAlertMessage, TelegramDataPackage

from ..alert.target_info import TargetInfo
from ..alert import AlertParsers
from ..config import CfgOption
from ..glade import GladeCatalog, GladeGalaxy
from ..json import JSON
from ..logger import log
from ..site import Site, Telescopes
from ..topic import TOPICS
from .. import site, sql

topics = CfgOption("topics", TOPICS, list)
timeout = CfgOption("timeout", -1, int)
start_date = CfgOption("start_date", datetime.now().isoformat(), datetime.fromisoformat)
output_dir = CfgOption("output_dir", "products", str)
default_site = CfgOption("default_site", "mondy_azt33ik", str)
ntasks = CfgOption("ntasks", 8, int)
max_localization_radius = CfgOption("max_localization_radius", 1, float)


if not os.path.exists(output_dir.get_value()):
    os.makedirs(output_dir.get_value())

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
    engine, session = sql.create_session()

    # Create alert table if it is not exist (e.g. fresh database)
    if not engine.has_table("alert"):
        sql.Alert.metadata.create_all(engine)

    # Unique hash for the alert message to not add duplicates to the database
    hash_md5 = md5(alert_msg.encode("utf-8")).hexdigest()

    # No errors, add the alert to the database
    loc = info.localization
    ra, dec = loc.center().ra.deg, loc.center().dec.deg
    r1 = r2 = loc.error_radius()
    localization = dumps(dict(r1=r1, r2=r2))
    alert_tab = sql.Alert(
        alert_message=alert_msg,
        ra_center=ra,
        dec_center=dec,
        error_radius1=r1,
        error_radius2=r2,
        localization=localization,
        trigger_date=info.trigger_date,
        event=info.event,
        origin=info.origin,
        importance=info.importance,
        md5=hash_md5,
    )

    with session:
        try:
            session.add(alert_tab)
        except Exception as e:
            pass

        try:
            session.commit()
        except Exception as e:
            log.warning("alert already in the database")
        else:
            log.info("alert added to the datebase")


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
    log.info("Reading messages from %s", start_date.get_value())
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
                "message has not passed validation, processing will not occur"
            )
            return

        que = self._data_queue

        text = message.value().decode("utf-8")
        log.info("\n" + "#" * 8 + "\n" + "New message\n" + str(datetime.now()))
        log.debug(text)

        # Retrieve target info
        log.debug("\nparsing alert ...")
        info = await parse_alert(message.value(), message.topic())
        if info is None:
            return
        log.debug(info)

        # # Save it to the local DB
        # log.debug("\ndumping alert info to database ...")
        # await dump_alert_to_db(text, info)

        # site = sites[default_site.get_value()]

        # Send the alert message before making plots and planning
        uuid = uuid4().hex

        async with self._lock:
            senter = info.origin
            body = (
                f"Event: {info.event}\n"
                f"Trigger: {info.trigger_date} UTC\n"
                f"Packet: {info.packet_type}\n"
                f"{info.localization.describe()}\n"
            )
            alert_message = TelegramAlertMessage(uuid, senter, body)
            await que.put(alert_message)

        sites = list(site.Telescopes[name] for name in site.default_sites.get_value())
        now = Time(datetime.now().isoformat(), format="isot")
        loc = info.localization
        jsons: list[JSON] = []
        axes: list[Axes] = []
        for i, s in enumerate(sites):
            start, stop = s.nearest_observation_window(now)

            # Can not be observed due to its location
            if start is None or stop is None:
                continue

            # If nearest observational window is located after 1 day since
            # today, it is considered that targets are not observable
            if info.origin == "LVC" and start.jd - now.jd >= 1:
                sorted_targets = []
                comment = (
                    "Targets do not cross the horizon at "
                    f"{s.full_name} within 24 hours."
                )
            else:
                # TODO: Let localization itself decide if it is needed
                # to make observational list for each site or just one
                # universal list
                sorted_targets, comment = await loc.observe(s, start, stop)

            # Dump sorted targets to JSON payload
            if sorted_targets:
                for tgt in sorted_targets:
                    tgt.name = info.event
                jsons.append(JSON(sorted_targets, s, info, comment))

                # TODO: Let localization itself decide if it is needed
                # to make plots for each site or just one universal plot
                if info.origin == "LVC":
                    if not axes:
                        axes = [loc.plot(s, start, stop, sorted_targets)]
                else:
                    ax = loc.plot(s, start, stop, sorted_targets)
                    axes.append(ax)

            if info.origin == "LVC":
                break

        # If payload is not empty
        if jsons:
            if len(jsons) > 1:
                jsons[0].add_jsons(*jsons[1:])
            payload = jsons[0].to_string()
            fname = os.path.join(output_dir.get_value(), f"targets_{uuid}.json")
            await save_file(fname, payload)
        else:
            # Do not send empty payload
            fname = ""

        plot_fnames: list[str] = []
        for i, ax in enumerate(axes):
            plot_fname = os.path.join(
                output_dir.get_value(), f"targets_{uuid}_{i+1}.png"
            )
            fig = ax.get_figure()
            fig.savefig(plot_fname)
            plot_fnames.append(plot_fname)

        # vis.plot_objects_on_localization(
        #     objects=sorted_targets, localization=info.localization,
        #     output_filename=plot_fname
        # )

        async with self._lock:
            data_package = TelegramDataPackage(uuid, info, fname, plot_fnames)
            await que.put(data_package)

        self._commited_messages += 1
        if self._commited_messages % self._max_commit_messages:
            self._consumer.commit(asynchronous=True)

        # start, stop = await site.nearest_observation_window(
        #     Time(info.trigger_date)
        # )
        # log.info(start.isot)
        # log.info(stop.isot)
        # sorted_targets = await info.localization.observe(site, start, stop)
        # log.info(sorted_targets)
        # fname = ""
        # if sorted_targets:
        #     log.debug(
        #         "\n created observational order for %d of targets observable",
        #         len(sorted_targets)
        #     )
        #     log.debug("saving ordered list of galaxies to JSON-file")
        #     serializer = JSON_Serializer(sorted_targets, site, info)
        #     dumped_galaxies = serializer.serialize()

        #     fname = os.path.join(
        #         output_dir.get_value(), f"targets_{uuid}.json"
        #     )

        #     await save_file(fname, dumped_galaxies)

        # plot_fname = os.path.join(
        #     output_dir.get_value(), f"targets_{uuid}.png"
        # )
        # vis.plot_objects_on_localization(
        #     objects=sorted_targets,
        #     localization=info.localization,
        #     output_filename=plot_fname
        # )

        # que.put(DataPackage(info, fname, plot_fname), block=True)

    def to_thread(self) -> Thread:
        return Thread(target=self.run, daemon=True)
