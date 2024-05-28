"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
bot.py (c) 2024
Desc: description
Created:  2024-03-23
Modified: 2024-04-02
"""

import asyncio
import os
from typing import Any, Awaitable, Literal
from aiogram import Dispatcher
from aiogram.dispatcher.storage import BaseStorage
from aiogram.contrib.fsm_storage.files import JSONStorage, MemoryStorage
from aiogram_dialog import DialogRegistry
from aiomisc import Service
from aware.data import AlertMessage, DataPackage
from aware.telegram.callbacks import Callback
from aware.telegram.callbacks.startup import MessageLoopCallback
from aware.logger import log
from aware.telegram.handlers.setup import setup_handlers
from aware.telegram.dialogs import setup_dialogs
from aware.telegram.decorator import AWARETelegramBot
from aiohttp import BasicAuth


TOKEN = os.getenv("AWARE_TG_API_TOKEN")
PROXY_URL = os.getenv("AWARE_TG_PROXY_URL")


class TelegramBot(Service):
    def __init__(
        self,
        token: str = TOKEN,
        proxy_url: str = PROXY_URL,
        proxy_auth: BasicAuth | None = None,
        queue: asyncio.Queue[AlertMessage | DataPackage] = None,
        storage: BaseStorage = MemoryStorage(),
        on_startup: list[Callback] | Callback | None = None,
        on_shutdown: list[Callback] | Callback | None = None,
        skip_updates: bool = True,
        mode: Literal["polling", "webhook"] = "polling",
        mode_kwargs: dict | None = None,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.bot = AWARETelegramBot(token, proxy=proxy_url, proxy_auth=proxy_auth)
        self.dp = Dispatcher(self.bot, storage=storage)
        self.dr = DialogRegistry(self.dp)
        self.queue = queue or asyncio.Queue()
        self.on_startup = on_startup
        self.on_shutdown = on_shutdown
        self.skip_updates = skip_updates
        self.mode = mode
        self.mode_kwargs = mode_kwargs or {}

    async def _run_on_startup(
        self, on_startup: list[Awaitable] | Awaitable | None = None
    ):
        if on_startup is not None:
            if hasattr(on_startup, "__iter__"):
                for startup_cb in on_startup:
                    await startup_cb(self.dp)
            else:
                await on_startup(self.dp)
        else:
            await MessageLoopCallback(self.queue)(self.dp)

    async def start(self) -> Any:
        # Register all handlers
        setup_handlers(self.dp)

        # Register all dialogs
        setup_dialogs(self.dr)

        # Run callbacks on initialization
        await self._run_on_startup(self.on_startup)

        # Skip updates, so, bot will not sent messages after restoring connection
        if self.skip_updates:
            await self.dp.skip_updates()

        # Currently, only long polling is supported
        # TODO: support webhooks
        
        if self.mode == "polling":
            await self.dp.start_polling(**self.mode_kwargs)
        elif self.mode == "webhooks":
            raise NotImplementedError("mode %s currently not supported" % self.mode)
        else:
            raise ValueError("mode %s is not recognized" % self.mode)
