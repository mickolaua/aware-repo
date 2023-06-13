import atexit
from datetime import datetime
from time import time

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
def main(telegram: bool):
    start_time = time()
    log.info("AWARE has started at %s", datetime.now())
    app = Application(run_telegram=telegram)
    app.run()

    atexit.register(goodbye, start_time)


if __name__ == "__main__":
    main()
