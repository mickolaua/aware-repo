"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
notify_subs.py (c) 2024
Desc: notify subscribers about something
Created:  2024-08-13
Modified: !date!
"""

import asyncclick as click
from aware.logger import log
from aware.telegram.util import select_subscribers
import asyncio
from aware.telegram.bot import TOKEN
import aiohttp
import ssl

# bot = Bot(TOKEN)
API_TELEGRAM_URL = "https://api.telegram.org/bot{}/sendMessage"


async def send_message(
    chat_id: str, message: str, token: str
) -> aiohttp.ClientResponse:
    sslcontext = ssl.create_default_context(purpose=ssl.Purpose.CLIENT_AUTH)
    async with aiohttp.ClientSession() as session:
        async with session.post(
            url=API_TELEGRAM_URL.format(token),
            headers={"Content-Type": "application/json"},
            params={
                "chat_id": chat_id,
                "text": message,
            },
        ) as response:
            resp = await response.read()
            if response.ok:
                click.echo(f"Successfully sent message to {chat_id}")
            else:
                click.echo(
                    "Message was not send due to an error. See the server response: "
                    f"{resp}"
                )
            return resp


@click.command("Notify subscribers")
@click.option("-m", "--message", help="Message to send to subscribers")
@click.option("-l", "--list-subscribers", help="List subscribers", is_flag=True)
@click.option("-s", "--subscriber", help="Subscriber ids", multiple=True)
async def main(message: str, list_subscribers: bool, subscriber: list[str]):
    subs = select_subscribers()
    if list_subscribers:
        click.echo(
            f"{'chat_id':^10} {'alert_type':^10} {'content_type':^10} "
            f"{'telescopes':^10}"
        )
        click.echo(
            f"{'-'*len('chat_id'):^10} {'-'*len('alert_type'):^10} "
            f"{'-'*len('content_type'):^10} {'-'*len('telescopes'):^10}"
        )
        for subscriber in subs:
            alerts = ",".join(subscriber["alert_type"])
            contents = ",".join(subscriber["content_type"])
            scopes = ",".join(subscriber["telescopes"])
            click.echo(
                f'{str(subscriber["chat_id"]):^10} '
                f"{alerts:^10} "
                f"{contents:^10} "
                f"{scopes:^10}"
            )
    else:
        if message:
            sub_ids = set([sub["chat_id"] for sub in subs])
            whitelist = set([int(s) for s in subscriber]) or sub_ids
            for sub_id in whitelist:
                if sub_id in sub_ids:
                    log.info("Sending message to subscriber with chat id %d", sub_id)
                    await send_message(sub_id, message, TOKEN)
                else:
                    log.warning(
                        "Subscriber with chat id %d is absent in the database", sub_id
                    )
        else:
            log.warning("Message is not provided")


if __name__ == "__main__":
    asyncio.run(main())
