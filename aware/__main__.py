import asyncio
import atexit
from datetime import datetime
from time import time

from .consumer import prepare_consumer, datetime_to_offset, ConsumeLoop
from .credentials import get_credentials
from .logger import log
from .telegram import create_tg_thread
from .app import Application
from .config import dev, cfg

 
def goodbye(start_time: float):
    log.info("\n--" + "\nGood Bye, See You Again!")
    log.info("I worked for %.3g min %.3g s.", *divmod(time() - start_time, 60))


def main():
    start_time = time()
    log.info("AWARE has started at %s", datetime.now())
    if dev.get_value():
        config = {
            'group.id': 'lvk', 
            # Infinitely read messages in development mode since the beginning of the 
            # partition
            'auto.offset.reset': 'earliest',
            'enable.auto.commit': False,
            'enable.auto.offset.store': True
        }
    else:
        config = {
            'group.id': 'lvk', 
            'enable.auto.commit': False,

            # Set maximum idle time out between reading topic records to 1 day (in msec)
            # to prevent application maximum poll interval exceeded
            'max.poll.interval.ms': 86_400_000
        }

    credits = get_credentials()
    consumer = prepare_consumer(config, credits)
    que = asyncio.Queue()
    consume_loop = ConsumeLoop(consumer, que)
    alert_thread = consume_loop.to_thread()
    tg_thread = create_tg_thread(que, skip_updates=True)

    app = Application([tg_thread, alert_thread])
    app.run()

    atexit.register(goodbye, start_time)


if __name__ == "__main__":
    main()

