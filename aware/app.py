from __future__ import annotations

import asyncio
from datetime import datetime
from typing import Literal

import aiomisc
import matplotlib as mpl
import uvloop

from aware.service import ServiceFactory

from .config import dev
from .sql.models import create_session
from .sql.util import create_alert_tables


asyncio.set_event_loop_policy(uvloop.EventLoopPolicy())


class Application:
    def __init__(
        self,
        name: str = "aware",
        run_telegram: bool = False,
        run_socket_server: bool = False,
        no_plot_gui=True,
        mode: Literal["prod", "test"] = "test",
        threads: int = 1,
    ):
        self._que = asyncio.Queue()
        self.name = name
        self.start_date = datetime.now()
        self._old_backend = mpl.get_backend()
        self.run_telegram = run_telegram
        self.run_socket_server = run_socket_server
        self.mode = mode

        if no_plot_gui:
            self._disable_gui_plot()

        self.loop = aiomisc.new_event_loop()
        self.threads = threads

    def _disable_gui_plot(self):
        mpl.use("agg")

    def restore_gui_plot(self):
        mpl.use(self._old_backend)

    def curr_plot_backend(self) -> str:
        return mpl.get_backend()

    def _shutdown(self):
        for task in asyncio.all_tasks(self.loop):
            try:
                task.cancel()
            except:
                ...
        self.loop.close()

    def _create_services(self):
        self.services = []
        factory = ServiceFactory(self._que)

        if self.mode == "test":
            self.services.append(factory.create_test_consumer())
        elif self.mode == "prod":
            self.services.append(factory.create_prod_consumer())

        if self.run_telegram:
            self.services.append(factory.create_telegram_bot())

        if self.run_socket_server:
            self.services.append(factory.create_socket_server())

    def run(self) -> float:
        engine, _ = create_session()
        create_alert_tables(engine)

        self._create_services()
        with aiomisc.entrypoint(*self.services, pool_size=self.threads) as loop:
            loop.run_forever()

        # Elapsed time
        time_delta = datetime.now() - self.start_date

        return time_delta.total_seconds()
