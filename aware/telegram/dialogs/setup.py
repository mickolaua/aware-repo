"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
setup.py (c) 2024
Desc: setup dialogs
Created:  2024-03-24
Modified: 2024-03-24
"""


from aiogram_dialog import DialogRegistry
from aware.telegram import dialogs


def setup_dialogs(registry: DialogRegistry):
    registry.register(dialogs.subscription_dialog)
    