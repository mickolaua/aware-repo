"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
unsubscribe.py (c) 2024
Desc: handler for unsubscribing from the bot
Created:  2024-03-24
Modified: 2024-03-24
"""

import aiogram

from aware.cache import pop_cache
from aware.sql.models import Subscriber, dbconnect
from aware.telegram.util import select_subscribers


async def unsubscribe_handler(m: aiogram.types.Message):
    for sub in select_subscribers():
        if sub["chat_id"] == m.chat.id:
            with dbconnect() as conn:
                obj = (
                    conn.query(Subscriber)
                    .filter(Subscriber.chat_id == m.chat.id)
                    .first()
                )
                conn.delete(obj)
            
            pop_cache("subscribers")
            return await m.reply("You was unsubscribed!")

    return await m.reply("You already unsubscribed!")