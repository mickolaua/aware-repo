from __future__ import annotations

import asyncio
from datetime import datetime
from queue import Queue
from threading import Lock, Thread
from typing import Any, Sequence

from .config import dev
from .consumer.main import ConsumeLoop, prepare_consumer
from .credentials import get_credentials
from .sql.models import create_session
from .sql.util import create_alert_tables
from .telegram.main import create_tg_thread
import matplotlib as mpl


CONSUMER_CONFIG = {
    "group.id": "lvk" if not dev.value else "lvk-test",
    # Infinitely read messages in development mode since the beginning of
    # the partition
    "auto.offset.reset": "earliest" if dev.get_value() else "latest",
    "enable.auto.commit": True,
    # Set maximum idle time out between reading topic records to 1 day
    # (in msec) to prevent application maximum poll interval exceeded
    "max.poll.interval.ms": 86_400_000,
}


class Application:
    def __init__(
        self,
        name: str = "aware",
        run_telegram: bool = True,
        no_plot_gui=True,
        consumer_config: dict[str, Any] | None = None,
    ):
        self._que = asyncio.Queue()
        self.name = name
        self.start_date = datetime.now()
        self._threads: list[Thread] = []
        self.run_telegram = run_telegram
        self._old_backend = mpl.get_backend()
        self.consumer_config = CONSUMER_CONFIG
        if consumer_config is not None:
            self.consumer_config.update(consumer_config)

        if no_plot_gui:
            self._disable_gui_plot()

    def _disable_gui_plot(self):
        mpl.use("agg")

    def restore_gui_plot(self):
        mpl.use(self._old_backend)

    def curr_plot_backend(self) -> str:
        return mpl.get_backend()

    def _spawn_consumer_thread(self):
        credits = get_credentials()
        consumer = prepare_consumer(CONSUMER_CONFIG, credits)
        consume_loop = ConsumeLoop(consumer, self._que)
        self._threads.append(consume_loop.to_thread())

    def _spawn_telegram_thread(self):
        main_tg_thread = create_tg_thread(que=self._que)
        self._threads.append(main_tg_thread)

    def _run_app_threads(self):
        # Split start and join of threads in different loops to not block
        # other threads from running
        for t in self._threads:
            t.start()

        for t in self._threads:
            t.join()

    def _spawn_threads(self):
        self._spawn_consumer_thread()
        if self.run_telegram:
            self._spawn_telegram_thread()

    def run(self) -> float:
        engine, _ = create_session()
        create_alert_tables(engine)
        self._spawn_threads()
        self._run_app_threads()

        # Elapsed time
        time_delta = datetime.now() - self.start_date

        return time_delta.total_seconds()
