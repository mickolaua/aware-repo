from __future__ import annotations

from queue import Queue
from threading import Lock, Thread
from typing import Sequence


class Application:
    def __init__(self, threads: Sequence[Thread]):
        self.threads = list(threads)
        self._lock = Lock()
        self._queue = Queue()
        
    def run(self) -> None:
        # Split start and join of threads in different loops to not block 
        # other threads from running
        for t in self.threads:
            t.start()

        for t in self.threads:
            t.join()
