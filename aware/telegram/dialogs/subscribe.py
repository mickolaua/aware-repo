"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
subscribe.py (c) 2024
Desc: dialog for the subscription to the bot
Created:  2024-03-24
Modified: 2024-03-24
"""

import asyncio
import operator
from threading import Thread
from typing import Any

import aiogram
from aiogram.dispatcher.filters.state import State, StatesGroup
from aiogram_dialog import (
    ChatEvent,
    Dialog,
    DialogManager,
    DialogRegistry,
    StartMode,
    Window,
)
from aiogram_dialog.widgets.kbd import (
    Back,
    Button,
    Cancel,
    Checkbox,
    ManagedCheckboxAdapter,
    Multiselect,
    Next,
    Row,
    Select,
    SwitchTo,
    Column,
)
from aiogram_dialog.widgets.text import Const, Format
from aiogram_dialog.widgets.when import Whenable


from aware import consumer, site
from aware.config import CfgOption, dev
from aware.logger import log
from aware.telegram.util import add_subscriber
from aware.visualization.find_chart import plot_find_chart
from aware.topic import (
    full_topic_name_to_short,
    short_topic_name_to_full,
    get_topics,
    GCNTOpic,
)
from aware.cache import pop_cache, write_cache, read_cache


class Form(StatesGroup):
    content_types = State()
    topics = State()
    sites = State()
    sub = State()


async def get_content_types(**kwargs):
    return {"content_types": [("Alerts", 1), ("Schedules", 2)]}


content_types_kbd = Multiselect(
    Format("✓ {item[0]}"),
    Format("{item[0]}"),
    id="s_content_types",
    item_id_getter=operator.itemgetter(
        0
    ),  # each item is a tuple with id on a first position
    items="content_types",  # we will use items from window data at a key `fruits`
)


async def get_data(**kwargs):
    topics = (
        (GCNTOpic.get_topic(name), i + 1)
        for i, name in enumerate(consumer.main.topics.value)
    )
    topic_names = []
    for topic, i in topics:
        if topic is not None:
            topic_names.append(topic.short_name)

    if topic_names:
        topic_names.sort()
        count = len(topic_names)
    else:
        topic_names = []
        count = 0

    return {
        "topics": zip(topic_names, range(len(topic_names))),
        "count": count,
    }


async def on_topic_selected(
    c: aiogram.types.CallbackQuery, widget: Any, manager: DialogManager, item_id: str
):
    log.debug("Topic selected: %s", item_id)
    # log.debug("Data: %s", manager.data["widget_data"]["s_content_types"])


topics_kbd = Multiselect(
    Format("✓ {item[0]}"),
    Format("{item[0]}"),
    id="s_topics",
    item_id_getter=operator.itemgetter(
        0
    ),  # each item is a tuple with id on a first position
    items="topics",  # we will use items from window data at a key `topics`
    on_click=on_topic_selected,
)


async def on_site_selected(
    c: aiogram.types.CallbackQuery, widget: Any, manager: DialogManager, item_id: str
):
    log.debug("Site selected: %s", item_id)
    # log.debug("Data: %s", manager.data["widget_data"]["s_sites"])


async def get_sites(**kwargs):
    sites = (
        (site.Telescopes.get(t_name, None), i + 1)
        for i, t_name in enumerate(site.default_sites.value)
    )
    site_names = []
    for tel, i in sites:
        if tel is not None:
            site_names.append(tel.full_name)

    if site_names:
        site_names.sort()
        count = len(site_names)
    else:
        site_names = []
        count = 0

    return {
        "sites": zip(site_names, range(len(site_names))),
        "count": count,
    }


sites_kbd = Multiselect(
    Format("✓ {item[0]}"),
    Format("{item[0]}"),
    id="s_sites",
    item_id_getter=operator.itemgetter(
        0
    ),  # each item is a tuple with id on a first position
    items="sites",  # we will use items from window data at a key `fruits`
    on_click=on_site_selected,
)


async def on_subscription(
    c: aiogram.types.CallbackQuery,
    widget: Any,
    manager: DialogManager,
):
    data = manager.data["aiogd_context"].widget_data
    await manager.done()

    log.debug("Data before finishing subscription: %s", data)
    content_types = data["s_content_types"]
    topic_by_short_name = {t.short_name: t for t in GCNTOpic.list_topics()}
    alert_topics = [
        topic_by_short_name[short_name].name for short_name in data["s_topics"]
    ]

    sites = [
        t_name
        for i, t_name in enumerate(site.default_sites.value)
        if site.Telescopes[t_name].full_name in data["s_sites"]
    ]

    ok = await add_subscriber(
        str(c.message.chat.id), content_types, alert_topics, sites
    )
    if ok:
        await asyncio.sleep(2)
        await c.message.reply("You are subscribed!")
        pop_cache("subscribers")
    else:
        await asyncio.sleep(2)
        await c.message.reply("You are NOT subscribed.\n")


subscription_dialog = Dialog(
    Window(
        Format("Choose what to receive"),
        content_types_kbd,
        Next(),
        Cancel(),
        state=Form.content_types,
        getter=get_content_types,
    ),
    Window(
        Format("Choose interesting topics"),
        Column(topics_kbd),
        Row(Back(), Next()),
        Cancel(),
        state=Form.topics,
        getter=get_data,
    ),
    Window(
        Format("Choose for which telescopes you want receive plots and schedules"),
        Column(sites_kbd),
        Row(Back(), Next()),
        Cancel(),
        state=Form.sites,
        getter=get_sites,
    ),
    Window(
        Format("Press `Subscribe` to finish subscribtion process"),
        Row(Back(), Button(Const("Subscribe"), id="b_sub", on_click=on_subscription)),
        Cancel(),
        state=Form.sub,
    ),
)
