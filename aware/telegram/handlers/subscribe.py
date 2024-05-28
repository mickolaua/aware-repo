"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
subscribe.py (c) 2024
Desc: description
Created:  2024-03-24
Modified: 2024-03-24
"""

import aiogram
from aiogram_dialog import DialogManager, StartMode

from aware.logger import log
from aware.telegram.util import select_subscribers
from aware.telegram.dialogs import subscription_dialog


async def subscribe_handler(m: aiogram.types.Message, dialog_manager: DialogManager):
    # Important: always set `mode=StartMode.RESET_STACK` you don't want to stack dialogs
    log.debug("User is trying to subscribe to the bot")

    ok = True
    for sub in select_subscribers():
        try:
            if sub["chat_id"] == m.chat.id:
                await m.reply("You already subscribed!")
                ok = False
        except Exception as e:
            ok = False
            log.error("cannot check if the user subscribed", exc_info=e)
            await m.reply("Sorry, cannot subsribed you due to internal error.")

    if ok:
        for state in subscription_dialog.states:
            if state.state.endswith("content_types"):
                return await dialog_manager.start(
                    state, mode=StartMode.RESET_STACK
                )
