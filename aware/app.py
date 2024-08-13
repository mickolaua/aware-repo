from __future__ import annotations

import asyncio
from datetime import datetime
from typing import Literal

import aiomisc
import matplotlib as mpl
import uvloop

from aware.logger import log
from aware.service import ServiceFactory
from aiomisc import Service

from .config import dev
from .sql.models import create_session
from .sql.util import create_alert_tables
from aware.data import AlertMessage, DataPackage
from collections import deque


class QueueHub(Service):
    """
    QueueHub is performing distribution of items from a main queue between
    queues owned by other clients. It is useful because when a multiple clients use
    the same queue, the one which gets an item first, actually removes it from the
    queue.
    """

    def __init__(
        self,
        master_queue: asyncio.Queue[AlertMessage | DataPackage],
        queues: list[asyncio.Queue[AlertMessage | DataPackage]] = None,
        **kwargs,
    ):
        self.master_queue = master_queue
        self.queues = deque(queues or [])

    async def distribute(self):
        """
        Perform actual distribution of items between master queue and other queues

        Raises
        ------
        ValueError
            if the queue is empty
        """
        if self.queues is None:
            raise ValueError(
                "No queues to distribute; use add_queue method to add queues"
            )

        while self.queues:
            try:
                item = await self.master_queue.get()
            except asyncio.QueueEmpty as e:
                log.error(
                    "Error occurred while sending item to hub memebers", exc_info=e
                )
            else:
                for q in self.queues:
                    await q.put(item)

    async def start(self):
        """Start the hub."""
        await self.distribute()

    async def add_queue(self, q: asyncio.Queue):
        """Add a queue to the hub

        Parameters
        ----------
        q : asyncio.Queue
            a queue to add
        """
        async with asyncio.Lock():
            self.queues.append(q)

    async def cleanup(self):
        """
        Cleanup the hub
        """
        async with asyncio.Lock():
            self.queues.clear()


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

        # if self.mode == "test":
        #     self.services.append(factory.create_test_consumer())
        # elif self.mode == "prod":
        #     self.services.append(factory.create_prod_consumer())

        # if self.run_telegram:
        #     self.services.append(factory.create_telegram_bot())

        # if self.run_socket_server:
        #     self.services.append(factory.create_socket_server())

        queues = []

        # Consumer fills the master queue with data
        if self.mode == "test":
            self.services.append(factory.create_test_consumer(queue=self._que))
        elif self.mode == "prod":
            self.services.append(factory.create_prod_consumer(queue=self._que))

        if self.run_telegram:
            tg_que = asyncio.Queue()
            self.services.append(factory.create_telegram_bot(queue=tg_que))
            queues.append(tg_que)

        if self.run_socket_server:
            socket_que = asyncio.Queue()
            self.services.append(factory.create_socket_server(queue=socket_que))
            queues.append(socket_que)

        hub = QueueHub(self._que, queues)
        self.services.insert(0, hub)

    def run(self) -> float:
        engine, _ = create_session()
        create_alert_tables(engine)

        self._create_services()
        with aiomisc.entrypoint(*self.services, pool_size=self.threads) as loop:
            loop.run_forever()

        # Elapsed time
        time_delta = datetime.now() - self.start_date

        return time_delta.total_seconds()
