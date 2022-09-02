import atexit
from datetime import datetime
from time import time
from .consumer import consume_loop
from gcn_kafka import Consumer
from .config import CLIENT_ID, CLIENT_SECRET
from .config import cfg
from .logger import log
import os


def goodbye(start_time: float) -> None:
    log.info("\n--")
    log.info("Good Bye, See You Again!")
    log.info("I worked for %.3g min %.3g s.", *divmod(time() - start_time, 60))


def main():
    start_time = time()
    log.info(
        "[%s] AWARE has started at %s", os.getlogin(), datetime.now()
    )
    config = {'group.id': 'read_last_unread',
              'auto.offset.reset': 'earliest',
              'enable.auto.commit': False}
    consumer = Consumer(config=config, 
                        client_id=CLIENT_ID, 
                        client_secret=CLIENT_SECRET)

    consume_loop(
        consumer, 
        cfg["consumer"]["topics"], 
        start_date=datetime.fromisoformat(cfg["consumer"]["start_date"])
    )

    atexit.register(goodbye, start_time)


if __name__ == "__main__":
    main()

