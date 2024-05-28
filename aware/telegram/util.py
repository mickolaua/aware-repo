from __future__ import annotations

import asyncio
from decimal import Decimal
import decimal
from hashlib import md5

from aiogram import Bot, types

from aware.cache import read_cache, write_cache
from aware.sql.models import Subscriber, dbconnect

from ..logger import log

# All messages with a size over this are considered long by Telegram
MAX_LENGTH_SHORT_MESSAGE = 512

# The Telegram flood control limits maximal count of messages per second (30 msg / sec)
MAX_MESSAGES_PER_MIN_GROUP = 20


def user_info_string(message: types.Message) -> str:
    usr_info = (
        f"{message.chat.username}@{message.chat.user_url} in chat {message.chat.id}"
    )
    return usr_info


async def send_short_messages_from_long_message(
    bot: Bot,
    long_message: types.Message | str,
    chat_id: str,
    max_len_short_msg=MAX_LENGTH_SHORT_MESSAGE,
    timeout_secs: float = 0.0,
    **kwargs,
):
    # Get bytes representation of the long message to determine its size
    if isinstance(long_message, str):
        bytes_msg = long_message.encode("utf-8")

    elif isinstance(long_message, types.Message):
        bytes_msg = long_message.text.encode("utf-8")

    else:
        raise ValueError(
            f"long_message must be of type {types.Message.__name__} "
            f"or {str.__name__}"
        )

    # Size of message in bytes
    size_msg = len(bytes_msg)

    # Determine the number of short messages from splitting the long message
    # in chunks of max_len_short_msg bytes or shorter
    n_msg, res = divmod(size_msg, max_len_short_msg)
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
        else (
            1.1 * MAX_MESSAGES_PER_MIN_GROUP / 60
            if n_msg > MAX_MESSAGES_PER_MIN_GROUP
            else 0.0
        )
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
        s.commit()

    return True


def select_subscribers():
    # key = "subscribers"
    # subscribers = read_cache(key)
    # if subscribers:
    #     subscribers = [Subscriber(**kws) for kws in subscribers]
    # else:
    sub_kwargs = []
    with dbconnect() as conn:
        subscribers = conn.query(Subscriber).yield_per(100)
        for i, sub in enumerate(subscribers):
            sub_kwargs.append(
                {
                    "chat_id": sub.chat_id,
                    "alert_type": sub.alert_type.split(","),
                    "content_type": sub.content_type.split(","),
                    "telescopes": sub.telescopes.split(","),
                }
            )
        # write_cache(key, sub_kwargs)

    # yield from subscribers
    return sub_kwargs


def _create_md5_file_checksum(filename: str) -> str:
    with open(filename, "rb") as f:
        file_data = f.read()
    file_hash = md5(file_data, usedforsecurity=False).hexdigest()
    return file_hash


def read_file_id(filename: str):
    file_hash = _create_md5_file_checksum(filename)
    key = f"file_{file_hash}"
    file_id = read_cache(key)
    return file_id


def write_file_id(filename, file_id):
    file_hash = _create_md5_file_checksum(filename)
    key = f"file_{file_hash}"
    write_cache(key, file_id, expire=3 * 86400)


async def add_subscriber(
    id: str, content_types: list[str], alert_topics: list[str], sites: list[str]
) -> bool:
    status = False

    def error_handler(exc_value: BaseException):
        log.error("user cannot subscribe: %s", exc_value, exc_info=exc_value)

    with dbconnect(exception_handler=error_handler) as s:
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
        status = True

    return status
