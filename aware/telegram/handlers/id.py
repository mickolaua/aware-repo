"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
id.py (c) 2024
Desc: description
Created:  2024-03-24
Modified: 2024-03-24
"""
import aiogram


async def user_id_handler(message: aiogram.types.Message):
    await message.reply(f"Your id is {message.chat.id}.")