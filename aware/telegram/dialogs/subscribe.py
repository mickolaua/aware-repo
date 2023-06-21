from __future__ import annotations

import asyncio
import operator
from typing import Any

import aiogram
from aiogram.dispatcher import FSMContext
from aiogram.dispatcher.filters import Text
from aiogram.dispatcher.filters.state import State, StatesGroup
from aiogram.utils import exceptions
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
)
from aiogram_dialog.widgets.text import Const, Format
from aiogram_dialog.widgets.when import Whenable

from ...logger import log
from .. import consumer, site


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
    topics = [
        (t.lstrip("gcn.classic.voevent."), i + 1)
        for i, t in enumerate(consumer.topics.get_value())
    ]
    return {
        "topics": topics,
        "count": len(topics),
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
    items="topics",  # we will use items from window data at a key `fruits`
    on_click=on_topic_selected,
)


async def on_site_selected(
    c: aiogram.types.CallbackQuery, widget: Any, manager: DialogManager, item_id: str
):
    log.debug("Site selected: %s", item_id)
    # log.debug("Data: %s", manager.data["widget_data"]["s_sites"])


async def get_sites(**kwargs):
    sites = [(t_name, i + 1) for i, (t_name, t) in enumerate(site.Telescopes.items())]
    return {
        "sites": sites,
        "count": len(sites),
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

    log.debug("Data before finish: %s", data)
    content_types = data["s_content_types"]
    alert_topics = ["gcn.classic.voevent." + t for t in data["s_topics"]]
    sites = data["s_sites"]
    ok = await add_subscriber(c.message.chat.id, content_types, alert_topics, sites)
    if ok:
        await asyncio.sleep(2)
        await c.bot.send_message(c.message.chat.id, "You are subscribed!")
    else:
        await asyncio.sleep(2)
        await c.bot.send_message(
            c.message.chat.id, "You are NOT subscribed due to internal error."
        )
    # log.debug("Data: %s", manager.data["widget_data"]["s_sites"])


dialog = Dialog(
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
        topics_kbd,
        Row(Back(), Next()),
        Cancel(),
        state=Form.topics,
        getter=get_data,
    ),
    Window(
        Format("Choose for which telescopes you want receive plots and schedules"),
        sites_kbd,
        Row(Back(), Next()),
        Cancel(),
        state=Form.sites,
        getter=get_sites,
    ),
    Window(
        Format("Almost done"),
        Row(Back(), Button(Const("Subscribe"), id="b_sub", on_click=on_subscription)),
        Cancel(),
        state=Form.sub,
    ),
)
