# from __future__ import annotations

# from .main import dp
# from aiogram import types
# from ..consumer import topics
# from ..site import default_sites
# from sqlalchemy import insert
# from ..sql.models import Subscriber
# from aiogram.dispatcher import FSMContext
# from aiogram.dispatcher.filters.state import State, StatesGroup


# class Form(StatesGroup):
#     content_types = State()
#     topics = State()
#     sites = State()
#     sub = State()
    

# @dp.message_handler(commands=['subscribe'])
# async def subscribe_user(message: types.Message):
#     """
#     Conversation's entry point
#     """
#     can_be_sent = False
#     try:
#         admins = await message.chat.get_administrators()
#         for admin in admins:
#             if message.from_user.id == admin.user.id:
#                 can_be_sent = True
#                 break
        
#     if not can_be_sent:
#         await message.reply("You have no rights to subscribe bot in this chat")
#         return

#     # Set state
#     await Form.content_types.set()

#     kb = types.InlineKeyboardMarkup()
#     kb.add(types.InlineKeyboardButton("Alerts"))
#     kb.add(types.InlineKeyboardButton("Schedules"))

#     await message.reply("What content you want to receive?", reply_markup=kb)


# @dp.message_handler(lambda message: not message.text.isdigit(), state=Form.content_types)
# async def incorrect_content_type(message: types.Message, state: FSMContext):
#     await message.reply("You did not choose any content type")
#     await state.finish()


# @dp.message_handler(state=Form.content_types)
# async def process_content_type(message: types.Message, state: FSMContext):
#     async with state.proxy() as data:
#         data['content_types'] = message.text.split("\n")

#         if "Alerts" in data["content_types"]:
#             await Form.topics.set()
#             kb = types.InlineKeyboardMarkup()
#             for t in topics.get_value():
#                 kb.add(types.InlineKeyboardButton(t))
#             return await message.reply("What topics are interesting for you?", reply_markup=kb)

#         elif "Schedules" in data["content_types"]:
#             await Form.sites.set()
#             kb = types.InlineKeyboardMarkup()
#             for t in default_sites.get_value():
#                 kb.add(types.InlineKeyboardButton(t))
#             return await message.reply("What telescopes are interesting for you?", reply_markup=kb)
        

# @dp.message_handler(lambda message: not message.text.isdigit(), state=Form.topics)
# async def incorrect_topics(message: types.Message, state: FSMContext):
#     await message.reply("You did not choose any topics")
#     await state.finish()


# @dp.message_handler(state=Form.topics)
# async def process_topics(message: types.Message, state: FSMContext):
#     async with state.proxy() as data:
#         data['alert_topics'] = message.text.split("\n")

#         if "Schedules" in data["content_types"]:
#             await Form.sites.set()
#             kb = types.InlineKeyboardMarkup()
#             for t in default_sites.get_value():
#                 kb.add(types.InlineKeyboardButton(t))
#             return await message.reply("What telescopes are interesting for you?", reply_markup=kb)
#         else:
#             await Form.sub.set()
#             await message.reply(f"You subscribed! Your settings bellow:")


# @dp.message_handler(lambda message: not message.text.isdigit(), state=Form.topics)
# async def incorrect_sites(message: types.Message, state: FSMContext):
#     await message.reply("You did not choose any telescopes")
#     await state.finish()


# @dp.message_handler(state=Form.sites)
# async def process_sites(message: types.Message, state: FSMContext):
#     async with state.proxy() as data:
#         data['sites'] = message.text.split("\n")

#         await Form.sub.set()
#         await message.reply(f"You subscribed! Your settings bellow:")


# @dp.message_handler(state=[Form.sub])
# async def process_sub(message: types.Message, state: FSMContext):
#     await state.finish()

#     async with state.proxy() as data:
#         if data["alert_topics"]:
#             settings_message = "\n[Topics]\n" + "\n".join(data["alert_topics"])
#         else:
#             settings_message = ""

#         if data["sites"]:
#             settings_message += "\n[Telescopes]\n" + "\n".join(data["sites"])

#         await message.reply(settings_message)

        