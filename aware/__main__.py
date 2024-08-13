import atexit
from datetime import datetime
import os
import sys
from time import time

from aware.config import save_cfg

from .logger import log
from .app import Application
import click


def goodbye(start_time: float):
    log.info("\n--" + "\nGood Bye, See You Again!")
    log.info("I worked for %.3g min %.3g s.", *divmod(time() - start_time, 60))


@click.command()
@click.option(
    "-t",
    "--telegram",
    is_flag=True,
    default=False,
    help="Run the Telegram thread",
)
@click.option(
    "-s",
    "--sockserver",
    is_flag=True,
    default=False,
    help="Run the socket server thread",
)
@click.option(
    "-m",
    "--mode",
    default="test",
    help="Application mode",
    type=click.Choice(["test", "prod"]),
)
@click.option("-T", "--threads", default=1, help="Number of threads in pool")
@click.option(
    "-d",
    "--dump",
    is_flag=True,
    help="Dump configuration to file (set by AWARE_CONFIG_FILE "
    "environment variable)",
)
def main(telegram: bool, sockserver: bool, mode: str, threads: int, dump: str):
    # if sockserver and telegram:
    #     log.error(
    #         "Simultaneous execution of Telegram bot and socket server not supported"
    #     )
    #     return 1

    if dump:
        if not os.getenv("AWARE_CONFIG_FILE", ""):
            log.error(
                "set environmental variable AWARE_CONFIG_FILE to the config file "
                "location"
            )
            return 1
        from aware import load_all_modules
        load_all_modules()
        save_cfg()
    else:
        start_time = time()
        log.info("AWARE has started at %s", datetime.now())
        app = Application(
            run_telegram=telegram,
            run_socket_server=sockserver,
            mode=mode,
            threads=threads,
        )
        app.run()
        atexit.register(goodbye, start_time)

    return 0


if __name__ == "__main__":
    sys.exit(main())
