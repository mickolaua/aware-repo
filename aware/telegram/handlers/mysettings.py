"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
mysettings.py (c) 2024
Desc: description
Created:  2024-03-24
Modified: 2024-03-24
"""

import aiogram
from aware import site
from aware.logger import log
from aware.telegram.util import select_subscribers
from aware.topic import GCNTOpic


async def mysettings_handler(message: aiogram.types.Message):
    log.debug("user looking for subsription settings")
    for sub in select_subscribers():
        if sub["chat_id"] == message.chat.id:
            reply = ""
            if "Alerts" in sub["content_type"] and sub["alert_type"]:
                reply += "[Alert types]\n"
                reply += "\n".join(
                   GCNTOpic.get_topic(t).short_name for t in sub["alert_type"]
                )

            if "Schedules" in sub["content_type"]:
                reply += "\n\n[Telescopes]\n"
                reply += "\n".join(
                    site.Telescopes[name].full_name for name in sub["telescopes"]
                )

            return await message.reply(reply)

    return await message.reply("You are not subscribed!")