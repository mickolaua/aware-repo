from __future__ import annotations

import asyncio
import os
import re
from datetime import datetime
from io import BytesIO
from multiprocessing import Process
from threading import Thread
from typing import Any

import aiogram
import numpy as np
import pandas as pd
from aiogram.contrib.fsm_storage.files import JSONStorage
from aiogram.dispatcher import FSMContext
from aiogram.dispatcher.filters import Text
from aiogram.dispatcher.filters.state import State, StatesGroup
from aiogram.utils import exceptions
from astroplan import FixedTarget
from astroplan.plots.finder import plot_finder_image
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.skyview import SkyView
from matplotlib import pyplot as plt

from aware.data import TelegramAlertMessage, TelegramDataPackage

from .. import consumer, site, sql
from ..config import CfgOption, dev
from ..logger import log

# Options
send_delay = CfgOption("send_delay", 2.0, float)

TOKEN = str(os.getenv("AWARE_TG_API_TOKEN", ""))
PROXY_URL = int(os.getenv("AWARE_TG_PROXY_URL", 0))
PROXY_AUTH = str(os.getenv("AWARE_TG_PROXY_AUTH", ""))

# All messages with a size over this are considered long by Telegram
# the flood control will count them and could initiate retry error
# (typically 40 sec)
MAX_LENGTH_SHORT_MESSAGE = 512
MAX_MESSAGES_PER_MIN_GROUP = 20


storage = JSONStorage("storage.json")
bot = aiogram.Bot(token=TOKEN, connections_limit=9999)
dp = aiogram.Dispatcher(bot, storage=storage)
queue = asyncio.Queue()
lock = asyncio.Lock()
running = False


def user_info_string(message: aiogram.types.Message) -> str:
    usr_info = (
        f"User {message.chat.username} @{message.chat.user_url} in chat "
        f"{message.chat.id}"
    )
    return usr_info


async def send_short_messages_from_long_message(
    bot: aiogram.Bot,
    long_message: aiogram.types.Message | str,
    chat_id: str,
    max_len_short_msg=MAX_LENGTH_SHORT_MESSAGE,
    timeout_secs: float = 0.0,
    **kwargs,
):
    # Get bytes representation of the long message to determine its size
    if isinstance(long_message, str):
        bytes_msg = long_message.encode("utf-8")

    elif isinstance(long_message, aiogram.types.Message):
        bytes_msg = long_message.text.encode("utf-8")

    else:
        raise ValueError(
            f"long_message must be of type {aiogram.types.Message.__name__} "
            f"or {str.__name__}"
        )

    # Size of message in bytes
    size_msg = len(bytes_msg)

    # Determine the number of short messages from splitting the long message
    # in chunks of max_len_short_msg bytes or shorter
    n_msg, res = divmod(size_msg, max_len_short_msg)
    log.debug("N = %d; R = %d", n_msg, res)
    if n_msg > 1:
        if res:
            n_msg += 1
    else:
        n_msg = 1

    # We should not hit the limit of 30 messages per minute or RetryError will occur.
    # So sleep before sending each message if needed.
    timeout_secs = (
        timeout_secs
        if timeout_secs
        else 1.1 * MAX_MESSAGES_PER_MIN_GROUP / 60
        if n_msg > MAX_MESSAGES_PER_MIN_GROUP
        else 0.0
    )

    # Sent short messages
    for i in range(n_msg):
        start = i * MAX_LENGTH_SHORT_MESSAGE
        stop = start + MAX_LENGTH_SHORT_MESSAGE + 1
        to_send = bytes_msg[start:stop].decode("utf-8")
        log.debug("run(): sending message %s", to_send)

        if timeout_secs:
            await asyncio.sleep(timeout_secs)

        await bot.send_message(chat_id, to_send, **kwargs)


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
    await message.reply("Wait for list of telescopes created")

    try:
        msg = "\n".join(
            [
                site.Telescopes[site_id].full_name
                for site_id in site.default_sites.get_value()
            ]
        )
    except Exception as e:
        msg = ""

    if not msg:
        await message.bot.send_message(
            message.chat.id, "Can not retrieve list of available telescopes"
        )
    else:
        # await message.reply(msg)
        try:
            await send_short_messages_from_long_message(
                message.bot, msg, message.chat.id
            )
        except:
            pass


@dp.message_handler(commands="status")
async def status(message: aiogram.types.Message):
    log.debug(
        "telegram.status(): user asked for status in chat (id=%d) on %s",
        message.chat.id,
        datetime.now(),
    )
    await message.reply(f"I am {'running' if running else 'idle'}")


@dp.message_handler(commands="run")
async def run(message: aiogram.types.Message):
    log.debug(
        "telegram.run(): ran in chat ID=%d on %s", message.chat.id, datetime.now()
    )
    global queue
    global lock
    global running

    async with lock:
        running = True

    await processing_loop(message, queue, lock)


@dp.message_handler(commands="stop")
async def stop(message: aiogram.types.Message):
    log.debug(
        "telegram.run(): requested to stop processing in chat ID=%d on %s",
        message.chat.id,
        datetime.now(),
    )
    global lock
    global running
    loop = asyncio.get_running_loop()
    loop.create_future()

    async with lock:
        if not running:
            return await message.reply("I am already idle")
        running = False

    await message.reply("I has stopped")


@dp.message_handler(regexp="[\d\w]+[\t ]\d+(\.\d+)?[\t ]\d+(\.\d+)?")
async def finding_chart(message: aiogram.types.Message):
    log.debug(
        "telegram.finding_chart(): user asked for finding chart in chat %s",
        message.chat.id,
    )
    await message.reply("Please, wait while DSS finding chart created")

    try:
        name, ra, dec = "", None, None
        for m in re.finditer("^[\w\d]+", message.text):
            name = m.string
            break

        ra_dec_match = [
            m.group()
            for m in re.finditer(r"(\d+)(\.\d+)?", message.text)
            if m is not None
        ]
        if ra_dec_match and len(ra_dec_match) == 2:
            ra, dec = ra_dec_match
            ra = float(ra)
            dec = float(dec)
    except Exception as e:
        log.error(
            "telegram.finding_chart(): failed to parse target name and coordinates "
            "from user message: %s",
            message.text,
            exc_info=e,
        )

    if name is None or ra is None or dec is None:
        return await message.reply(
            "Can not figure out how to extract target name and its coordinate from "
            "your message."
        )

    try:
        target = FixedTarget(SkyCoord(ra, dec, unit=["deg"] * 2), name=name)
        img = SkyView.get_images(
            target.coord,
            survey="DSS",
            coordinates="J2000",
            radius=15 * u.arcmin,
            scaling="Log",
        )[0]
        hdr = img[0].header
        data = img[0].data
        vmin = data.mean() - data.std()
        vmax = data.mean() + data.std()
        wcs = WCS(hdr)
        plt.subplot(projection=wcs)
        plt.imshow(data, vmin=vmin, vmax=vmax, cmap="gray_r")
        plt.plot(
            target.coord.ra,
            target.coord.dec,
            ls="",
            marker="*",
            color="red",
            label=name,
            transform=plt.gca().get_transform("world"),
        )
        plt.xlabel("RA")
        plt.ylabel("Dec")
        plt.legend()
        io = BytesIO()
        plt.savefig(io)
        io.seek(0)

        await message.bot.send_photo(
            message.chat.id,
            aiogram.types.InputFile(io),
            reply_to_message_id=message.message_id,
        )
    except Exception as e:
        log.error(
            "telegram.finding_chart(): can not create finding chart for target %s "
            "at (%.5f, %.5f)",
            ra,
            dec,
            exc_info=e,
        )
        await message.reply(
            "Sorry, unable to create the finding chart with provided parameters"
        )


@dp.message_handler(commands="help")
async def help(message: aiogram.types.Message):
    log.debug(
        "telegram.help(): user asked for help in chat ID=%d on %s",
        message.chat.id,
        datetime.now(),
    )
    await message.reply(
        """
/help - Ask for help
/run - Start alert receiving and observation planning
/stop - Stop this process
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
        "telegram.run(): requested list of topics in chat ID=%d on %s",
        message.chat.id,
        datetime.now(),
    )
    await message.reply(
        "Please, wait for the list available topics: (could be sent in "
        "many messages )"
    )
    response_message = "List of available topics:" + "\n".join(
        consumer.topics.get_value()
    )
    await send_short_messages_from_long_message(
        message.bot, response_message, message.chat.id
    )


# @functools.lru_cache(maxsize=8)
asyncio.coroutine


async def processing_loop(
    message: aiogram.types.Message,
    queue: asyncio.Queue[TelegramDataPackage | TelegramAlertMessage],
    lock: asyncio.Lock,
):
    await message.reply("I'm starting broadcasting alerts and planning observations.")
    global running
    while running:
        # Wait until data came then unpack it
        async with lock:
            try:
                data = await queue.get()
                queue.task_done()
                log.debug(
                    "telegram.run(): received item %s created on %s with id %s",
                    data,
                    data.created,
                    data.id,
                )
            except Exception as e:
                log.error("telegram.run(): failed to fetch data from queue", exc_info=e)
                continue

        if isinstance(data, TelegramAlertMessage):
            try:
                msg = (
                    f"From: {data.senter}\n"
                    f"{data.body}\n"
                    "Wait for the plots and observational plans in upcoming messages "
                    "(if any)."
                )
                await send_short_messages_from_long_message(
                    message.bot, msg, message.chat.id
                )
            except Exception as e:
                log.error("telegram.run(): can not send messages", exc_info=e)
            continue

        elif isinstance(data, TelegramDataPackage):
            # Sent available plots (e.g. visibility plots or localization plots)
            info = data.target_info
            if data.plot_fnames:
                msg = f"\nPlot filename(s) for {info.origin} {info.event} below:\n"
                msg += "\n".join([os.path.basename(fn) for fn in data.plot_fnames])

                # Split the message into short ones and sent them
                await asyncio.sleep(1)

                try:
                    await send_short_messages_from_long_message(
                        message.bot, msg, message.chat.id
                    )
                except Exception as e:
                    log.error("telegram.run(): can not send messages", exc_info=e)

                # Send localizaion map plots
                for fn in data.plot_fnames:
                    try:
                        await message.bot.send_photo(
                            message.chat.id, aiogram.types.InputFile(fn)
                        )
                    except Exception as e:
                        log.error("telegram.run(): can not send photo", exc_info=e)

            # If available send observational targets sorted in optimal order
            if data.json_filename:
                msg = (
                    f"Observational plan for {info.origin} {info.event} is in the file "
                    f"{os.path.basename(data.json_filename)})"
                )
                await asyncio.sleep(1)
                try:
                    await message.bot.send_message(message.chat.id, msg)
                except Exception as e:
                    log.error("telegram.run(): can not send messages", exc_info=e)

                await asyncio.sleep(1)

                try:
                    await message.bot.send_document(
                        message.chat.id, aiogram.types.InputFile(data.json_filename)
                    )
                except Exception as e:
                    log.error(
                        "telegram.run(): can not send JSON file of targets", exc_info=e
                    )
            continue


def start_polling(
    dp: aiogram.Dispatcher,
    que: asyncio.Queue[TelegramDataPackage],
    **kwargs: dict[Any, Any],
):
    # Windows-related stuff
    WINDOWS = os.name.startswith("win")
    if WINDOWS:
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())
        kwargs["loop"] = asyncio.new_event_loop()

    # Run event loop inside worker or polling would never find the event loop
    loop = kwargs.pop("loop", asyncio.new_event_loop())
    asyncio.set_event_loop(loop)
    global queue
    queue = que

    # Start long polling; note, loop already removed from kwargs
    aiogram.executor.start_polling(dp, loop=loop, **{k: v for k, v in kwargs.items()})


def create_tg_thread(*args, **kwargs) -> Thread:
    return Thread(
        target=start_polling,
        args=(dp, *args),
        name="tg_bot_thread",
        kwargs=kwargs,
        daemon=True,
    )
