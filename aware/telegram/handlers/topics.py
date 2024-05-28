"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
topics.py (c) 2024
Desc: List all available GCN topics 
Created:  2024-03-24
Modified: 2024-03-24
"""
import aiogram
from aware import consumer
from aware.logger import log
from aware.telegram.util import send_short_messages_from_long_message
from aware.topic import GCNTOpic


async def topics_handler(message: aiogram.types.Message):
    log.debug(
        "requested list of topics",
    )
    response_message = "Available topics:\n" + "\n".join(
        GCNTOpic.get_topic(t).short_name for t in consumer.topics.get_value()
    )
    await send_short_messages_from_long_message(
        message.bot, response_message, message.chat.id
    )