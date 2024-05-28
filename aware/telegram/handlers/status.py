"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
status.py (c) 2024
Desc: Show bot status
Created:  2024-03-24
Modified: 2024-03-24
"""
import aiogram

from aware.logger import log


async def status_handler(message: aiogram.types.Message):
    log.debug("user asked for status")
    await message.reply("I am serving.")