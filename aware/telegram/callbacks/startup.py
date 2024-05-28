import asyncio
import os
import string
from functools import partial
from textwrap import dedent

import aiogram
from aiogram import Dispatcher
from aiogram.types import ParseMode
from aiogram.utils import markdown
from emoji import emojize

from aware.cache import read_cache
from aware.data import AlertMessage, DataPackage
from aware.logger import log
from aware.telegram.callbacks.callback import Callback
from aware.telegram.util import read_file_id, select_subscribers, write_file_id
from aware.util import render_number_as_rich_utf8
from aware.topic import full_topic_name_to_short, GCNTOpic

# Monkey-patch function, to avoid provding alias each time
emojize = partial(emojize, language="alias")


def alert_message_to_rich_text(alert_message: AlertMessage) -> str:
    """
    Create a rich text alert message to be sent by the Telegram bot
    """
    info = alert_message.target_info
    loc = info.localization
    date = info.trigger_date.date()
    time = info.trigger_date.time().replace(microsecond=0)
    raw_tag = f"#{info.origin.replace(' ', '')}{info.event.replace(' ', '')}"
    tag_letters = []
    allowed_letters = set(string.ascii_letters)
    allowed_letters.update(string.digits)
    for letter in raw_tag:
        if letter in allowed_letters:
            tag_letters.append(letter)

    tag = "".join(tag_letters)

    locdescr = info.describe() if info.localization else ""
    metadescr = info.description
    first_mention = (
        alert_message.first_mention
        if alert_message.first_mention
        else f"{alert_message.senter} {info.event}"
    )
    msg = emojize(
        markdown.text(
            ":new:" if alert_message.is_new else "",
            ":boom: :email: Received alert message from GCN NASA\n",
            "\n:satellite:",
            markdown.escape_md(f" Experiment: {alert_message.senter}\n"),
            ":id:",
            markdown.escape_md(f" Trigger ID or event name: {info.event}\n"),
            ":date:",
            markdown.escape_md(f" Trigger date and time: {date} {time} UT\n"),
            ":1st_place_medal:",
            markdown.escape_md(f" First mention: {first_mention}\n"),
            ":label:",
            markdown.escape_md(f" Packet: {alert_message.alert_type}\n"),
            "\n:point_down:" if locdescr else "",
            markdown.escape_md(f" Localization parameters below:") if locdescr else "",
            ":point_down:\n\n" if locdescr else "",
            markdown.escape_md(locdescr) if locdescr else "",
            "\n:arrows_counterclockwise:" if alert_message.update else "",
            (
                markdown.escape_md(
                    f" Localization for {alert_message.first_mention} refined\n"
                )
                if alert_message.update and loc
                else ""
            ),
            "\n:point_down:" if metadescr else "",
            markdown.escape_md(f" Event metadata below:") if metadescr else "",
            ":point_down:\n\n" if metadescr else "",
            markdown.escape_md(metadescr) if metadescr else "",
            "\n\n:x: EVENT HAS BEEN CANCELLED :x:\n" if info.rejected else "",
            (
                ":exclamation: DO NOT OBSERVE THIS AND CROSS\-MATCHED EVENTS :exclamation:\n"
                if info.rejected
                else ""
            ),
            "\n\n:pushpin:" if tag else "",
            markdown.escape_md(f" Tag: #{tag}") if tag else "",
            sep="",
        )
    )

    return msg


def data_package_to_rich_plot_text(data_package: DataPackage) -> str:
    """Create a rich formatted plot message to be sent by the Telegram bot"""
    msg = emojize(
        markdown.text(
            ":earth_africa: Observational plot\n",
            ":telescope:",
            markdown.escape_md(f" Site: {data_package.site.full_name}\n"),
            ":satellite:",
            markdown.escape_md(f" Origin: {data_package.target_info.origin}\n"),
            ":id:",
            markdown.escape_md(f" Event: {data_package.target_info.event}\n"),
            ":page_facing_up:",
            markdown.escape_md(
                f"Filename: {os.path.basename(data_package.plot_fname)}"
            ),
        )
    )
    return msg


def data_package_to_rich_plan_text(data_package: DataPackage) -> str:
    """Create a rich formatted observatio plan message to be sent by the Telegram bot"""
    msg = emojize(
        markdown.text(
            ":memo: Observational plan\n",
            ":telescope:",
            markdown.escape_md(f" Site: {data_package.site.full_name}\n"),
            ":satellite:",
            markdown.escape_md(f" Origin: {data_package.target_info.origin}\n"),
            ":id:",
            markdown.escape_md(f" Event: {data_package.target_info.event}\n"),
            ":page_facing_up:",
            markdown.escape_md(
                f"Filename: {os.path.basename(data_package.plan_filename)}"
            ),
            (
                "\n\n:speech_balloon: Comments to an observer :point_down:\n\n"
                if data_package.comment
                else ""
            ),
            markdown.escape_md(data_package.comment) if data_package.comment else "",
        )
    )
    return msg


def simple_alert_message(alert_message: AlertMessage) -> str:
    info = alert_message.target_info
    loc = info.localization
    date = info.trigger_date.date()
    time = info.trigger_date.time().replace(microsecond=0)
    raw_tag = f"#{info.origin.replace(' ', '')}{info.event.replace(' ', '')}"
    tag_letters = []
    allowed_letters = set(string.ascii_letters)
    allowed_letters.update(string.digits)
    for letter in raw_tag:
        if letter in allowed_letters:
            tag_letters.append(letter)

    tag = "".join(tag_letters)

    locdescr = info.describe().strip(" \n") if info.localization else ""
    metadescr = info.description.strip(" \n") if info.description else ""
    first_mention = (
        alert_message.first_mention or f"{alert_message.senter} {info.event}"
    )

    msg = "".join(
        [
            "GCN Alert\n--\n",
            f"Type: {GCNTOpic.get_topic(alert_message.alert_type).short_name}\n",
            f"Experiment: {alert_message.senter}\n",
            f"ID: {info.event}\n",
            f"First mention: {first_mention}\n",
            f"T0: {date} {time} UT\n",
            f"\nLocalization\n--\n{locdescr}" if locdescr else "",
            "\n" if locdescr else "",
            f"\nMeta\n--\n{metadescr}" if metadescr else "",
            "\n" if metadescr else "",
            (
                f"EVENT HAS BEEN CANCELLED, DO NOT OBSERVE THIS AND RELATED EVENTS!\n"
                if info.rejected
                else ""
            ),
            f"\n#{tag}",
        ]
    )
    return markdown.escape_md(msg)


def simple_plot_message(data_package: DataPackage) -> str:
    msg = dedent(
        f"""
        Observational plot
        --
        Event: {data_package.target_info.origin} {data_package.target_info.event}
        Telescope: {data_package.site.full_name}
        Filename: {os.path.basename(data_package.plot_fname)}
        """
    )
    return  markdown.escape_md(msg)


def simple_plan_message(data_package: DataPackage) -> str:
    msg = dedent(
        f"""
        Observational plan
        --
        Event: {data_package.target_info.origin} {data_package.target_info.event}
        Telescope: {data_package.site.full_name}
        Filename: {os.path.basename(data_package.plan_filename)}
        {"Comment:" + data_package.comment if data_package.comment else ""},
        """
    )
    return  markdown.escape_md(msg)


class MessageLoopCallback(Callback):
    def __init__(self, queue: asyncio.Queue[DataPackage | AlertMessage]):
        self.queue = queue

    async def __call__(self, dp: Dispatcher):
        """
        This is the main heavy task for the bot, sends alert messages, and observation
        plans to users. It utilizes the queue, which contains the necessary messages.

        Parameters
        ----------
        dp : Dispatcher
            a bot dispatcher (convinience for all aiogram startup callbacks)
        queue : asyncio.Queue[DataPackage  |  AlertMessage]
            a message queue (connected with planner)
        """
        log.debug("started broadcasting alerts and observational data")
        queue = self.queue

        # Since while True blocks other commands, we need to create a task
        async def inner():
            while True:
                # Wait until data came then unpack it
                try:
                    data = await queue.get()
                    queue.task_done()
                    log.debug(
                        "received item %s created on %s with id %s",
                        data,
                        data.created,
                        data.id,
                    )
                except Exception as e:
                    log.error("failed to fetch data from queue", exc_info=e)

                # Get bot
                bot = dp.bot

                # Send alert messages
                if isinstance(data, AlertMessage):
                    try:
                        for sub in select_subscribers():
                            if (
                                "Alerts" in sub["content_type"]
                                and data.alert_type in sub["alert_type"]
                            ):
                                # msg = alert_message_to_rich_text(data)
                                msg = simple_alert_message(data)
                                await asyncio.sleep(2)
                                await bot.send_message(
                                    sub["chat_id"],
                                    msg,
                                    parse_mode=ParseMode.MARKDOWN_V2,
                                )
                    except Exception as e:
                        log.error("can not send alert message: %s", e, exc_info=e)

                elif isinstance(data, DataPackage):
                    # Sent available plots (e.g. visibility plots or localization plots)
                    for sub in select_subscribers():
                        info = data.target_info
                        if (
                            "Schedules" in sub["content_type"]
                            and data.alert_type in sub["alert_type"]
                            and data.site.name in sub["telescopes"]
                        ):
                            if data.plot_fname:
                                # msg = emojize(
                                #     markdown.text(
                                #         ":map:",
                                #         "[Observational Plot]",
                                #         f"Site: {data.site.full_name}",
                                #         f"Origin: {info.origin}",
                                #         f"Event: {info.event},",
                                #         f"Filename: {os.path.basename(data.plot_fname)}",
                                #     )
                                # )
                                msg = simple_plot_message(data)

                                await asyncio.sleep(2)

                                try:
                                    # Check if file is already uploaded, then provide only
                                    # its id. This helps saving Telegram server resources
                                    # as well as our server recourses.
                                    plot_exists = False
                                    if read_cache(data.plot_fname):
                                        plot_file = read_file_id(data.plot_fname)
                                        plot_exists = True
                                    else:
                                        plot_file = aiogram.types.InputFile(
                                            data.plot_fname
                                        )

                                    plot_msg = await bot.send_photo(
                                        sub["chat_id"],
                                        photo=plot_file,
                                        caption=msg,
                                        parse_mode=ParseMode.MARKDOWN_V2,
                                    )
                                    if not plot_exists:
                                        write_file_id(
                                            data.plot_fname, plot_msg.photo[-1].file_id
                                        )
                                except Exception as e:
                                    log.error(
                                        "can not send observation plot: %s",
                                        e,
                                        exc_info=e,
                                    )

                            if data.plan_filename:
                                msg = simple_plan_message(data)

                                await asyncio.sleep(2)

                                try:
                                    # Send observational plan, but checking if it has 
                                    # sent already.
                                    plan_exists = False
                                    if read_cache(data.plan_filename):
                                        plan_file = read_file_id(data.plan_filename)
                                        plan_exists = True
                                    else:
                                        plan_file = aiogram.types.InputFile(
                                            data.plan_filename
                                        )

                                    plan_msg = await bot.send_document(
                                        sub["chat_id"],
                                        document=plan_file,
                                        caption=msg,
                                        parse_mode=ParseMode.MARKDOWN_V2,
                                    )
                                    if not plan_exists:
                                        write_file_id(
                                            data.plan_filename,
                                            plan_msg.document.file_id,
                                        )

                                except Exception as e:
                                    log.error(
                                        "can not send observation plan: %s",
                                        e,
                                        exc_info=e,
                                    )

        task = asyncio.ensure_future(inner())

        def done_callback(fut: asyncio.Future):
            try:
                fut.result()
            except:
                fut.cancel()

        task.add_done_callback(done_callback)
