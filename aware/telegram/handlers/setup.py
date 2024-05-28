"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
setup.py (c) 2024
Desc: Setup the bot handlers
Created:  2024-03-24
Modified: 2024-03-25
"""

import aiogram
from aiogram import Dispatcher

from aware.telegram import handlers


def setup_handlers(dp: Dispatcher):
    dp.register_message_handler(handlers.sites_handler, commands="telescopes")
    dp.register_message_handler(handlers.status_handler, commands="status")
    dp.register_errors_handler(
        handlers.retry_handler, exception=aiogram.exceptions.RetryAfter
    )
    dp.register_message_handler(
        handlers.visibility_handler,
        lambda msg: msg.text.startswith("/visplot") and len(msg.text.split(" ")) == 6,
        commands="visplot",
    )
    dp.register_message_handler(handlers.user_id_handler, commands="id")
    dp.register_message_handler(handlers.topics_handler, commands="topics")
    dp.register_message_handler(handlers.mysettings_handler, commands="mysettings")
    dp.register_message_handler(handlers.unsubscribe_handler, commands="unsub")
    dp.register_message_handler(handlers.subscribe_handler, commands="sub")
    dp.register_message_handler(handlers.visibility_handler, commands="visplot")
    dp.register_message_handler(
        handlers.findchart_handler,
        lambda message: message.text.startswith("/findchart")
        and len(message.text.split(" ")) == 3,
        commands="findchart",
    )
