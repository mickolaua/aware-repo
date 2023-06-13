import asyncio
import operator
import queue
from typing import Any

import aiogram
from aiogram.dispatcher.filters.state import State, StatesGroup
from aiogram.types import Message
from aiogram_dialog import Dialog, DialogManager, Window, DialogRegistry
from aiogram_dialog.tools import render_preview, render_transitions
from aiogram_dialog.widgets.input import MessageInput
from aiogram_dialog.widgets.kbd import (
    Back,
    Button,
    Cancel,
    Multiselect,
    Next,
    Row,
    Select,
    Column,
)
from aiogram_dialog.widgets.text import Const, Format
import pytest

from aware import consumer, site
from aware.logger import log
from aware.topic import full_topic_name_to_short, short_topic_name_to_full
from aware.telegram.main import bot, dp


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
        (full_topic_name_to_short(t), i + 1)
        for i, t in enumerate(consumer.topics.get_value())
    ]
    return {
        "topics": topics,
        "count": len(topics),
    }

async def on_topic_selected(
    c: aiogram.types.CallbackQuery, widget: Select, manager: DialogManager, item_id: str
):
    log.debug("Topic selected: %s", item_id)
    widget.text = Format(f"✓{widget.text}")
    log.debug(widget.text)
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
    item_id_getter=lambda l: l,  # each item is a tuple with id on a first position
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


dialog = Dialog(
    Window(
        Format("Choose what to receive"),
        content_types_kbd,
        Next(),
        Cancel(),
        state=Form.content_types,
        preview_data={"content_types": [("Alerts", 1), ("Schedules", 2)]},
    ),
    Window(
        Format("Choose interesting topics"),
        Column(topics_kbd),
        # Row(Back(), Next()),
        # Cancel(),
        state=Form.topics,
        preview_data={
            "topics": [
                (full_topic_name_to_short(t), i + 1)
                for i, t in enumerate(consumer.topics.get_value())
            ]
        },
    ),
    Window(
        Format("Choose for which telescopes you want receive plots and schedules"),
        Column(sites_kbd),
        Row(Back(), Next()),
        Cancel(),
        state=Form.sites,
        preview_data={
            "sites": [
                (site.Telescopes[t_name].full_name, i + 1)
                for i, t_name in enumerate(site.default_sites.value)
            ]
        },
    ),
    Window(
        Format("Almost done"),
        Row(Back(), Button(Const("Subscribe"), id="b_sub", on_click=on_subscription)),
        Cancel(),
        state=Form.sub,
    ),
)


@pytest.mark.skip
@pytest.mark.asyncio
async def test():
    reg = DialogRegistry(dp)
    reg.register(dialog)
    render_transitions([dialog])
    await render_preview(reg, "dialog.html", simulate_events=True)


if __name__ == "__main__":
    asyncio.run(test())
