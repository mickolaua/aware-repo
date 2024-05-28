"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
callback.py (c) 2024
Desc: description
Created:  2024-03-24
Modified: 2024-03-24
"""
from typing import Any

from aiogram import Dispatcher


class Callback:
    async def __call__(self, dp: Dispatcher):
        raise NotImplementedError
