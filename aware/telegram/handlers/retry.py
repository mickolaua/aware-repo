"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
retry.py (c) 2024
Desc: Handler for retry error
Created:  2024-03-24
Modified: 2024-03-24
"""
import asyncio

import aiogram


async def retry_handler(
    update: aiogram.types.Update, exception: aiogram.exceptions.RetryAfter
):
    await asyncio.sleep(exception.timeout + 3)