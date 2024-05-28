"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
list_telescopes.py (c) 2024
Desc: description
Created:  2024-03-24
Modified: 2024-03-24
"""

from datetime import datetime
import aiogram

from aware import site
from aware.cache import read_cache, write_cache
from aware.logger import log
from aware.telegram.util import send_short_messages_from_long_message


__all__ = ["show_available_sites"]


async def sites_handler(message: aiogram.types.Message):
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
