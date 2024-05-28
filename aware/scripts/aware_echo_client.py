from __future__ import annotations

import asyncio

import click as click
from aiomisc import entrypoint
from aiomisc.service import TCPClient


async def read_until_EOF(reader: asyncio.StreamReader, eof: bytes) -> str:
    data = await reader.readuntil(eof)
    return data.decode()


class SocketClient(TCPClient):

    async def handle_connection(
        self, reader: asyncio.StreamReader, writer: asyncio.StreamWriter
    ):
        """Handles a connection with a client. Note: `aiomisc.TCPClient` already handles 
        connection errors.

        Parameters
        ----------
        reader : asyncio.StreamReader
            _description_
        writer : asyncio.StreamWriter
            _description_
        """
        try:
            while True:
                msg = await read_until_EOF(reader, b"--STOP_PROGRAM--")
                click.echo(msg)
        finally:
            await writer.drain()
            writer.close()
            await writer.wait_closed()


def entry_point(host: str, port: int):
    with entrypoint(SocketClient(host, port)) as loop:
        loop.run_forever()


@click.command()
@click.argument("host")
@click.option("-p", "--port", type=int, default=8080, help="Default port")
def main(host: str, port: int):
    """
    Listen host:port for observation programs. The program displays all the
    messages, which have been received. Currently, no any post-processing of
    messages performed, but it can be implemented by end users.

    The message format (by rows):
    ```
    telescope_id
    metadata
    --START_PROGRAM--
    target1
    target2
    ...
    targetN
    --STOP_PROGRAM--
    ```
    """
    entry_point(host, port)


if __name__ == "__main__":
    main()
