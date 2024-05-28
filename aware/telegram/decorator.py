"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
decorator.py (c) 2024
Desc: Decorator for bot new features, like sending rich text messages
Created:  2024-03-24
Modified: 2024-03-24
"""

from aiogram import Bot
from aiogram import types
from aiogram.utils import markdown
from emoji import emojize


# Here we implement decoration using inheritance actually, because aiogram checks for
# type.
class AWARETelegramBot(Bot):
    """
    Decorator class for aiogram.Bot providing extended functionality.
    """

    # async def send_message(self, chat_id: int, text: str, **kwargs) -> types.Message:
    #     """Send message with rich formatting (markdown)"""

    #     try:
    #         del kwargs["parse_mode"]
    #     except KeyError:
    #         pass

    #     msg = markdown.text(emojize(text, language="alias"))
    #     return await super().send_message(
    #         chat_id,
    #         text=text,
    #         parse_mode=types.ParseMode.MARKDOWN_V2,
    #         **kwargs,
    #     )

    async def reply_error(self, text: str, message: types.Message) -> types.Message:
        return await message.reply(
            emojize(markdown.text(":x:", markdown.bold(text)), language="alias"),
            parse_mode=types.ParseMode.MARKDOWN_V2,
        )

    async def echo(self, message: types.Message) -> types.Message:
        return await self.send_message(
            message.chat.id, message.text, reply_to_message_id=message.id
        )
