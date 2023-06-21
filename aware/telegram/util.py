from __future__ import annotations

import asyncio

from aiogram import Bot, types

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
    engine, session = create_session()
    with session:
        stmt = select(Subscriber)
        subscribers: list[tuple[Subscriber]] = session.execute(stmt).fetchall()
        for i, (sub,) in enumerate(subscribers):
            sub.alert_type = sub.alert_type.split(",")
            sub.content_type = sub.content_type.split(",")
            sub.telescopes = sub.telescopes.split(",")
            subscribers[i] = sub

    # yield from subscribers
    return subscribers
