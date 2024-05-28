from __future__ import annotations

import asyncio
import operator
import os
import re
import warnings
from datetime import datetime
from io import BytesIO
from multiprocessing import Process
from threading import Thread
from typing import Any

import aiogram
import numpy as np
import pandas as pd
from aiogram.contrib.fsm_storage.files import JSONStorage
from aiogram.contrib.fsm_storage.redis import RedisStorage2
from aiogram.dispatcher import FSMContext
from aiogram.dispatcher.filters import Text
from aiogram.dispatcher.filters.state import State, StatesGroup
from aiogram.utils import exceptions
from aiogram_dialog import (
    ChatEvent,
    Dialog,
    DialogManager,
    DialogRegistry,
    StartMode,
    Window,
)
from aiogram_dialog.widgets.kbd import (
    Back,
    Button,
    Cancel,
    Checkbox,
    Column,
    ManagedCheckboxAdapter,
    Multiselect,
    Next,
    Row,
    Select,
    SwitchTo,
)
from aiogram_dialog.widgets.text import Const, Format
from aiogram_dialog.widgets.when import Whenable
from astroplan import FixedTarget
from astroplan.plots.finder import plot_finder_image
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.skyview import SkyView
from matplotlib import pyplot as plt
from sqlalchemy import exc, select

from aware.data import TelegramAlertMessage, TelegramDataPackage
from aware.sql.models import Subscriber, create_session

from .. import consumer, site, sql
from ..cache import pop_cache, read_cache, write_cache
from ..config import CfgOption, dev
from ..logger import log
from ..topic import full_topic_name_to_short, short_topic_name_to_full
from ..visualization.find_chart import plot_find_chart
from .util import send_short_messages_from_long_message, user_info_string

warnings.warn(
    "Using aware.telegram.main considered deprecated, please use "
    "aware.telegram.bot instead.",
    DeprecationWarning,
)


# Options
send_delay = CfgOption("send_delay", 2.0, float)

TOKEN = str(os.getenv("AWARE_TG_API_TOKEN", ""))
PROXY_URL = int(os.getenv("AWARE_TG_PROXY_URL", 0))
PROXY_AUTH = str(os.getenv("AWARE_TG_PROXY_AUTH", ""))


storage = JSONStorage("storage.json")
bot = aiogram.Bot(token=TOKEN, connections_limit=9999)
dp = aiogram.Dispatcher(bot, storage=storage)
registry = DialogRegistry(dp)
queue = asyncio.Queue()
lock = asyncio.Lock()
running = False


@dp.errors_handler(exception=aiogram.exceptions.RetryAfter)
async def retry_handler(
    update: aiogram.types.Update, exception: aiogram.exceptions.RetryAfter
):
    await asyncio.sleep(exception.timeout + 3)


@dp.message_handler(commands=["telescopes"], state="*")
async def show_available_sites(message: aiogram.types.Message):
    log.debug(
        "user requested list of available sites in chat ID=%d on %s",
        message.chat.id,
        datetime.now(),
    )
    key = "telescopes"
    try:
        telescopes = read_cache(key)
        if not telescopes:
            telescopes = sorted(
                [
                    site.Telescopes[site_id].full_name
                    for site_id in site.default_sites.get_value()
                ]
            )
            write_cache(key, telescopes)

        msg = "\n".join(telescopes)
    except Exception as e:
        msg = ""

    if not msg:
        await message.bot.send_message(
            message.chat.id, "Can not retrieve list of available telescopes"
        )
    else:
        try:
            await send_short_messages_from_long_message(
                message.bot, msg, message.chat.id
            )
        except:
            pass


@dp.message_handler(commands="status")
async def status(message: aiogram.types.Message):
    log.debug(
        "telegram.status(): user asked for status",
    )
    await message.reply(f"I am serving.")


@dp.message_handler(commands="visplot")
@dp.message_handler(
    lambda msg: msg.text.startswith("/visplot") and len(msg.text.split(" ")) == 6
)
async def status(message: aiogram.types.Message):
    log.debug(
        "telegram.status(): user asked for visibility plot",
    )
    try:
        cmd_name, name, ra, dec, date, telescope_id = message.text.split(" ")
    except TypeError:
        return await message.reply(
            "Invalid message format, should be like this: <name> <ra> <dec> <date> <telescope_id>"
        )

    if telescope_id not in site.Telescopes:
        await message.reply(f"Telescope with id={telescope_id} is not available")
    else:
        telescope = site.Telescopes[telescope_id]
        coord = SkyCoord(ra, dec)
        target = FixedTarget(coord, name)
        init_time = Time(date, format="isot")
        key = f"telescope_id={telescope_id} target={name} ra={coord.ra.deg:.0f} dec={coord.dec.deg:.0f} date={init_time.datetime.date()}"
        plot_file_content = read_cache(key)
        if plot_file_content is None:
            start_time, end_time = telescope.nearest_observation_window(init_time)
            ax = telescope.plot_airmass(
                [target], start_time=start_time, end_time=end_time
            )
            plot_file = BytesIO()
            ax.get_figure().savefig(plot_file)
            plot_file.seek(0)
            plot_file_content = plot_file.read()
            write_cache(key, plot_file_content)

        plot_file = BytesIO(plot_file_content)

        await message.bot.send_photo(
            message.chat.id,
            aiogram.types.InputFile(plot_file),
            caption="Your visibility plot",
            reply_to_message_id=message.message_id,
        )


# @dp.message_handler(
#     lambda message: message.text.startswith("/filter")
#     and len(message.text.split(" ")) == 2
# )
# async def set_filter(message: aiogram.types.Message):
#     _, msg_filter = message.text.split(" ")

#     can_filter = False
#     try:
#         admins = await message.chat.get_administrators()
#         ids = set([a.user.id for a in admins])
#         can_filter = message.from_id in ids and bot.id in ids
#     except aiogram.exceptions.BadRequest:
#         can_filter = True

#     if can_filter:
#         key = f"chat_{message.chat.id}_filters"
#         write_cache(key, msg_filter)


# @dp.message_handler()
# async def remove_message(message: aiogram.types.Message):
#     can_filter = False
#     try:
#         admins = await message.chat.get_administrators()
#         ids = set([a.user.id for a in admins])
#         can_filter = message.from_id in ids and bot.id in ids
#     except aiogram.exceptions.BadRequest:
#         can_filter = True

#     if can_filter:
#         msg_filter = read_cache(f"chat_{message.chat.id}_filters")
#         if (
#             msg_filter is not None
#             and message.from_id == message.bot.id
#             and re.search(msg_filter, message.text)
#         ):
#             await message.bot.delete_message(message.chat.id, message.message_id)


# @dp.message_handler(commands="run")
# async def run(message: aiogram.types.Message):
#     log.debug(
#         "telegram.run(): ran in chat ID=%d on %s", message.chat.id, datetime.now()
#     )
#     global queue
#     global lock
#     global running

#     async with lock:
#         running = True

#     await processing_loop(message, queue, lock)


# @dp.message_handler(commands="stop")
# async def stop(message: aiogram.types.Message):
#     log.debug(
#         "requested to stop processing in chat ID=%d on %s",
#         message.chat.id,
#         datetime.now(),
#     )

#     is_admin_or_creator = False
#     try:
#         admins = await message.chat.get_administrators()
#         for admin in admins:
#             if admin.status == aiogram.types.ChatMemberStatus.ADMINISTRATOR:
#                 is_admin_or_creator = True

#             if admin.status == aiogram.types.ChatMemberStatus.OWNER:
#                 is_admin_or_creator = True
#     except aiogram.exceptions.BadRequest:
#         # In private chat no admins
#         is_admin_or_creator = True

#     async def _shutdown():
#         global lock
#         global running

#         async with lock:
#             if not running:
#                 return await message.reply("I am already idle")
#             running = False

#         await message.reply("I has stopped")

#     if is_admin_or_creator:
#         asyncio.create_task(_shutdown())
#         await asyncio.sleep(0)
#     else:
#         await message.reply("You are not admin or owner of this chat!")


# @dp.message_handler(commands="resume")
# async def resume(message: aiogram.types.Message):
#     log.debug(
#         "requested to resume processing in chat ID=%d on %s",
#         message.chat.id,
#         datetime.now(),
#     )

#     is_admin_or_creator = False
#     try:
#         admins = await message.chat.get_administrators()
#         for admin in admins:
#             if admin.status == aiogram.types.ChatMemberStatus.ADMINISTRATOR:
#                 is_admin_or_creator = True

#             if admin.status == aiogram.types.ChatMemberStatus.OWNER:
#                 is_admin_or_creator = True
#     except aiogram.exceptions.BadRequest:
#         # In private chat no admins
#         is_admin_or_creator = True

#     async def _resume():
#         global lock
#         global running

#         async with lock:
#             if running:
#                 return await message.reply("I am already running")
#             running = True

#         await message.reply("I am resuming running")

#     if is_admin_or_creator:
#         asyncio.create_task(_resume())
#         await asyncio.sleep(0)
#     else:
#         await message.reply("You are not admin or owner of this chat!")


@dp.message_handler(commands=["findchart"])
@dp.message_handler(
    lambda message: message.text.startswith("/findchart")
    and len(message.text.split(" ")) == 3
)
async def finding_chart(message: aiogram.types.Message):
    log.debug(
        "user asked for finding chart",
    )
    try:
        _, name, ra, dec = message.text.split(" ")
        try:
            ra = float(ra)
            ra = Angle(ra, unit="deg")
        except ValueError:
            ra = Angle(ra, unit="hr")

        try:
            dec = float(dec)
        except ValueError:
            ...

        dec = Angle(dec, unit="deg")

        coord = SkyCoord(ra, dec)
        key = f"ra={coord.ra.deg:.0f}, dec={coord.dec.deg:.0f}, name={name}"
        plot_file_content = read_cache(key)
        if plot_file_content is None:
            ax = plot_find_chart(coord, 15 * u.arcmin, name)
            plot_file = BytesIO()
            ax.get_figure().savefig(plot_file, format="png")
            plot_file.seek(0)
            plot_file_content = plot_file.read()
            write_cache(key, plot_file_content)

        plot_file = BytesIO(plot_file_content)

        await message.bot.send_photo(
            message.chat.id,
            aiogram.types.InputFile(plot_file),
            reply_to_message_id=message.message_id,
        )
    except IndexError as e:
        log.error(
            "field is not available in SkyView or no connection: %s",
            e,
            exc_info=e,
        )
        await message.reply("Field is not available in Sky View")
    except Exception as e:
        log.error(
            "can not create finding chart, message was %s",
            message.text,
            exc_info=e,
        )
        await message.reply(
            "Unable to create the finding chart with provided parameters. Tip: check "
            "the format, it must be like so: /findchart <name> <ra> <dec> (brackets "
            "are ommited)."
        )


@dp.message_handler(commands=["id"])
async def get_ids(message: aiogram.types.Message):
    await asyncio.sleep(2)
    await message.reply(f"{message.chat.id}")


@dp.message_handler(commands="help")
async def help(message: aiogram.types.Message):
    log.debug(
        "user asked for help",
    )
    await message.reply(
        """
/help - Ask for help
/status - Display if processing is ran or not
/telescopes - Get list of supported telescopes
/topics - Get list of GCN events to be received

Text commands:
1. Finding chart image
    Type TRANSIENT_NAME 133.33 35.33 to get DSS finding chart of the target
"""
    )


@dp.message_handler(commands="topics")
async def topics(message: aiogram.types.Message):
    log.debug(
        "requested list of topics",
    )
    response_message = "Available topics:\n" + "\n".join(
        t.lstrip("gcn.classic.voevent.") for t in consumer.topics.get_value()
    )
    await send_short_messages_from_long_message(
        message.bot, response_message, message.chat.id
    )


async def add_subscriber(
    id: str, content_types: list[str], alert_topics: list[str], sites: list[str]
) -> bool:
    _, session = create_session()
    with session as s:
        sub = Subscriber()

        if alert_topics:
            if len(alert_topics) > 1:
                alert_topics = ",".join(alert_topics)
            else:
                alert_topics = alert_topics[0]
        else:
            alert_topics = ""

        if sites:
            if len(sites) > 1:
                sites = ",".join(sites)
            else:
                sites = sites[0]
        else:
            sites = ""

        if content_types:
            if len(sites) > 1:
                content_types = ",".join(content_types)
            else:
                content_types = content_types[0]
        else:
            content_types = ""

        sub.alert_type = alert_topics
        sub.telescopes = sites
        sub.content_type = content_types
        sub.chat_id = id

        s.add(sub)

        try:
            s.commit()
        except exc.IntegrityError as e:
            log.error("user cannot subscribe", exc_info=e)
            s.rollback()
            return False

    return True


def select_subscribers():
    key = "subscribers"
    subscribers = read_cache(key)
    if subscribers:
        subscribers = [Subscriber(**kws) for kws in subscribers]
    else:
        engine, session = create_session()
        sub_kwargs = []
        with session:
            subscribers = session.query(Subscriber).all()
            for i, sub in enumerate(subscribers):
                sub.alert_type = sub.alert_type.split(",")
                sub.content_type = sub.content_type.split(",")
                sub.telescopes = sub.telescopes.split(",")
                subscribers[i] = sub
                sub_kwargs.append(
                    {
                        "chat_id": sub.chat_id,
                        "alert_type": sub.alert_type,
                        "content_type": sub.content_type,
                        "telescopes": sub.telescopes,
                    }
                )
        write_cache(key, sub_kwargs)

    # yield from subscribers
    return subscribers


@dp.message_handler(commands=["mysettings"])
async def mysettings(message: aiogram.types.Message):
    for sub in select_subscribers():
        if sub.chat_id == message.chat.id:
            reply = ""
            if "Alerts" in sub.content_type and sub.alert_type:
                reply += "[Alert types]\n"
                reply += "\n".join(
                    t.lstrip("gcn.classic.voevent.") for t in sub.alert_type
                )

            if "Schedules" in sub.content_type:
                reply += "\n[Telescopes]\n"
                reply += "\n".join(
                    site.Telescopes[name].full_name for name in sub.telescopes
                )

            return await message.reply(reply)

    return await message.reply("You are not subscribed!")


@dp.message_handler(commands=["unsub"])
async def unsubscribe(m: aiogram.types.Message):
    for sub in select_subscribers():
        if sub.chat_id == m.chat.id:
            engine, session = create_session()
            with session:
                obj = (
                    session.query(Subscriber)
                    .filter(Subscriber.chat_id == m.chat.id)
                    .first()
                )
                session.delete(obj)
                session.commit()

            thread = Thread(
                None,
                pop_cache,
                args=("subscribers",),
                name="subscriber_unsub_cache_thread",
            )
            thread.start()
            thread.join()

            return await m.reply("You was unsubscribed!")

    return await m.reply("You already unsubscribed!")


@dp.message_handler(commands=["sub"])
async def subscribe(m: aiogram.types.Message, dialog_manager: DialogManager):
    # Important: always set `mode=StartMode.RESET_STACK` you don't want to stack dialogs
    log.debug("User is trying to subscribe to the bot")

    ok = True
    for sub in select_subscribers():
        try:
            if sub.chat_id == m.chat.id:
                await m.reply("You already subscribed!")
                ok = False
        except Exception as e:
            ok = False
            log.error("cannot check if the user subscribed", exc_info=e)
            await m.reply("Sorry, cannot subsribed you due to internal error.")

    if ok:
        await dialog_manager.start(Form.content_types, mode=StartMode.RESET_STACK)


class Form(StatesGroup):
    content_types = State()
    topics = State()
    sites = State()
    sub = State()


async def get_content_types(**kwargs):
    return {"content_types": [("Alerts", 1), ("Schedules", 2)]}


content_types_kbd = Multiselect(
    Format("✓ {item[0]}"),
    Format("{item[0]}"),
    id="s_content_types",
    item_id_getter=operator.itemgetter(
        0
    ),  # each item is a tuple with id on a first position
    items="content_types",  # we will use items from window data at a key `fruits`
)


async def get_data(**kwargs):
    topics = [
        (full_topic_name_to_short(t), i + 1)
        for i, t in enumerate(consumer.topics.get_value())
    ]
    return {
        "topics": topics,
        "count": len(topics),
    }


async def on_topic_selected(
    c: aiogram.types.CallbackQuery, widget: Any, manager: DialogManager, item_id: str
):
    log.debug("Topic selected: %s", item_id)
    # log.debug("Data: %s", manager.data["widget_data"]["s_content_types"])


topics_kbd = Multiselect(
    Format("✓ {item[0]}"),
    Format("{item[0]}"),
    id="s_topics",
    item_id_getter=operator.itemgetter(
        0
    ),  # each item is a tuple with id on a first position
    items="topics",  # we will use items from window data at a key `topics`
    on_click=on_topic_selected,
)


async def on_site_selected(
    c: aiogram.types.CallbackQuery, widget: Any, manager: DialogManager, item_id: str
):
    log.debug("Site selected: %s", item_id)
    # log.debug("Data: %s", manager.data["widget_data"]["s_sites"])


async def get_sites(**kwargs):
    sites = [
        (site.Telescopes[t_name].full_name, i + 1)
        for i, t_name in enumerate(site.default_sites.value)
    ]
    return {
        "sites": sites,
        "count": len(sites),
    }


sites_kbd = Multiselect(
    Format("✓ {item[0]}"),
    Format("{item[0]}"),
    id="s_sites",
    item_id_getter=operator.itemgetter(
        0
    ),  # each item is a tuple with id on a first position
    items="sites",  # we will use items from window data at a key `fruits`
    on_click=on_site_selected,
)


async def on_subscription(
    c: aiogram.types.CallbackQuery,
    widget: Any,
    manager: DialogManager,
):
    data = manager.data["aiogd_context"].widget_data
    await manager.done()

    log.debug("Data before finish: %s", data)
    content_types = data["s_content_types"]
    alert_topics = [short_topic_name_to_full(t) for t in data["s_topics"]]
    sites = [
        t_name
        for i, t_name in enumerate(site.default_sites.value)
        if site.Telescopes[t_name].full_name in data["s_sites"]
    ]

    ok = await add_subscriber(
        str(c.message.chat.id), content_types, alert_topics, sites
    )
    if ok:
        await asyncio.sleep(2)
        await c.message.reply("You are subscribed!")
        thread = Thread(
            None,
            pop_cache,
            args=("subscribers",),
            name="subscriber_expire_cache_thread",
        )
        thread.start()
        thread.join()
    else:
        await asyncio.sleep(2)
        await c.message.reply("You are NOT subscribed.\n")


dialog = Dialog(
    Window(
        Format("Choose what to receive"),
        content_types_kbd,
        Next(),
        Cancel(),
        state=Form.content_types,
        getter=get_content_types,
    ),
    Window(
        Format("Choose interesting topics"),
        Column(topics_kbd),
        Row(Back(), Next()),
        Cancel(),
        state=Form.topics,
        getter=get_data,
    ),
    Window(
        Format("Choose for which telescopes you want receive plots and schedules"),
        Column(sites_kbd),
        Row(Back(), Next()),
        Cancel(),
        state=Form.sites,
        getter=get_sites,
    ),
    Window(
        Format("Press `Subscribe` to finish subscribtion process"),
        Row(Back(), Button(Const("Subscribe"), id="b_sub", on_click=on_subscription)),
        Cancel(),
        state=Form.sub,
    ),
)

registry.register(dialog)


async def processing_loop(
    # message: aiogram.types.Message,
    _,
    queue: asyncio.Queue[TelegramDataPackage | TelegramAlertMessage] = queue,
    lock: asyncio.Lock = lock,
):
    # await message.reply("Started broadcasting alerts and observational data")
    log.debug("started broadcasting alerts and observational data")
    global running

    while running:
        # Wait until data came then unpack it
        async with lock:
            try:
                data = await queue.get()
                queue.task_done()
                log.debug(
                    "received item %s created on %s with id %s",
                    data,
                    data.created,
                    data.id,
                )
            except Exception as e:
                log.error("failed to fetch data from queue", exc_info=e)

        if isinstance(data, TelegramAlertMessage):
            try:
                for sub in select_subscribers():
                    if (
                        "Alerts" in sub.content_type
                        and data.alert_type in sub.alert_type
                    ):
                        msg = data.body
                        await asyncio.sleep(2)
                        await bot.send_message(sub.chat_id, msg)
                        # await send_short_messages_from_long_message(
                        #     message.bot, msg, message.chat.id
                        # )
            except Exception as e:
                log.error("can not send messages", exc_info=e)
        elif isinstance(data, TelegramDataPackage):
            # Sent available plots (e.g. visibility plots or localization plots)
            for sub in select_subscribers():
                info = data.target_info
                if (
                    "Schedules" in sub.content_type
                    and data.alert_type in sub.alert_type
                    and data.site.name in sub.telescopes
                ):
                    if data.plot_fname:
                        msg = (
                            "[Observational Plot]\n"
                            f"Site: {data.site.full_name}\n"
                            f"Origin: {info.origin}\n"
                            f"Event: {info.event}\n"
                            f"Filename: {os.path.basename(data.plot_fname)}"
                        )

                        # try:
                        #     await asyncio.sleep(2)
                        #     await bot.send_message(sub.chat_id, msg)
                        # except Exception as e:
                        #     log.error("can not send message", exc_info=e)

                        # Split the message into short ones and sent them
                        await asyncio.sleep(2)
                        try:
                            await asyncio.sleep(2)
                            await bot.send_photo(
                                sub.chat_id,
                                aiogram.types.InputFile(data.plot_fname),
                                caption=msg,
                            )
                        except Exception as e:
                            log.error("can not send photo", exc_info=e)

                    if data.plan_filename:
                        msg = (
                            f"[Observational plan]\n"
                            f"Site: {data.site.full_name}\n"
                            f"Origin: {info.origin}\n"
                            f"Event: {info.event}\n"
                            f"Filename: "
                            f"{os.path.basename(data.plan_filename)}\n"
                            f"{data.comment}\n"
                        )
                        # await asyncio.sleep(2)
                        # try:
                        #     await asyncio.sleep(2)
                        #     await bot.send_message(sub.chat_id, msg)
                        # except Exception as e:
                        #     log.error("can not send messages", exc_info=e)

                        await asyncio.sleep(2)

                        try:
                            await asyncio.sleep(2)
                            await bot.send_document(
                                sub.chat_id,
                                aiogram.types.InputFile(data.plan_filename),
                                caption=msg,
                            )
                        except Exception as e:
                            log.error(
                                "can not send file with targets",
                                exc_info=e,
                            )


def start_polling(
    dp: aiogram.Dispatcher,
    que: asyncio.Queue[TelegramAlertMessage | TelegramDataPackage],
    **kwargs: dict[Any, Any],
):
    # Windows-related stuff
    WINDOWS = os.name.startswith("win")
    if WINDOWS:
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())
        kwargs["loop"] = asyncio.new_event_loop()

    # Run event loop inside worker or polling would never find the event loop
    loop: asyncio.AbstractEventLoop = kwargs.pop("loop", asyncio.new_event_loop())
    asyncio.set_event_loop(loop)
    global queue
    queue = que

    # Start long polling; note, loop already removed from kwargs
    async def on_startup(dp, queue=queue, lock=lock):
        global running
        running = True
        asyncio.create_task(processing_loop(dp, queue, lock))
        await asyncio.sleep(0)

    aiogram.executor.start_polling(
        dp,
        loop=loop,
        **{k: v for k, v in kwargs.items()},
        skip_updates=True,
        on_startup=on_startup,
    )


def create_tg_thread(*args, **kwargs) -> Thread:
    return Thread(
        target=start_polling,
        args=(dp, *args),
        name="telegram",
        kwargs=kwargs,
        daemon=True,
    )
