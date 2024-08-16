from __future__ import annotations

import asyncio

import click as click
from aiomisc import entrypoint
from aiomisc.service import TCPClient

async def readline(reader: asyncio.StreamReader) -> str:
    data = await reader.readuntil()
    return data.decode()


class SocketClient(TCPClient):

    async def handle_connection(
        self, reader: asyncio.StreamReader, writer: asyncio.StreamWriter
    ):
        try:
            while True:
                msg = await readline(reader)
                click.echo(msg)
        finally:
            await writer.drain()
            writer.close()
            await writer.wait_closed()


def entry_point(host: str, port: int):
    with entrypoint(SocketClient(host, port)) as loop:
        try:
            loop.run_forever()
        except KeyboardInterrupt:
            exit()


@click.command()
@click.argument("host", default="127.0.0.1")
@click.option("-p", "--port", type=int, default=55555, help="Default port")
def main(host: str, port: int):
    """
    Listens host:port for alert messages and observation programs.
    """
    entry_point(host, port)


if __name__ == "__main__":
    main()
