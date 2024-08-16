import asyncio
import os

import aiofiles
from aiomisc.service import TCPServer
import numpy as np

from aware.config import CfgOption
from aware.data import AlertMessage, DataPackage
from aware.logger import log
import re


IPV4_PORT_REGEX = (
    # First 3 numbers of IP: make sure each number <= 255
    r"^((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}"
    # Last 3 numbers of IP
    r"((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?))"
    # Match port below 60000
    r"(:([0-9]|[1-9][0-9]|[1-9][0-9]{2}|[1-9][0-9]{3}|[1-5][0-9][0-9]{3}|"
    # Match port between 60000 and 65535
    r"[6][0-5][0-5][0-3][0-5]))?$"
)


hostname = CfgOption("hostname", "127.0.0.1", str, comment="Hostname of the server")
port = CfgOption("port", 55555, int, comment="Port of the server")
max_connections = CfgOption(
    "max_connections",
    5,
    lambda x: max(1, x),
    comment="Maximum number of connections to the server",
)
send_alert_message = CfgOption(
    "send_alert_message",
    False,
    bool,
    comment="Send alert message via socket connection?",
)
send_cancelled_alerts = CfgOption(
    "send_cancelled_alerts",
    False,
    bool,
    comment="Send messages on cancelled alerts via socket connection?",
)
client_name_filters = CfgOption(
    "client_name_filters",
    [IPV4_PORT_REGEX],
    lambda x: [re.compile(i) for i in x],
    comment="Regular expression filters to validate client ip:port against",
)


PROGRAM_HEADER = "--START_PROGRAM--"
PROGRAM_FOOTER = "--STOP_PROGRAM--"


async def form_message(data: AlertMessage | DataPackage) -> bytes:
    if isinstance(data, DataPackage):
        if data.target_info.rejected and not send_cancelled_alerts.value:
            return b""
        else:
            if os.path.exists(data.plan_filename) and os.path.isfile(
                data.plan_filename
            ):
                lines = []
                async with aiofiles.open(data.plan_filename, "rb") as f:
                    async for line in f:
                        lines.append(line.rstrip().decode())

                header = f"{data.site.name}\n{PROGRAM_HEADER}\n"
                body = "\n".join(lines)
                footer = f"\n{PROGRAM_FOOTER}\n"

                return (header + body + footer).encode("ascii")
            else:
                return b""
    elif isinstance(data, AlertMessage):
        if not send_alert_message.value or (
            data.target_info.rejected and not send_cancelled_alerts.value
        ):
            return b""
        else:
            return data.body.encode("ascii")
    else:
        return b""


def is_allowed(ip: str, port: int) -> bool:
    """Check if the given ip and port are allowed to connect to the server

    Parameters
    ----------
    ip : str
        an IP address
    port : int
        a port number

    Returns
    -------
    bool
        True if the given ip and port are allowed
    """
    allowed = False
    for i, f in enumerate(client_name_filters.value):
        if re.fullmatch(f, f"{ip}:{port}"):
            allowed = True
            break

    return allowed


async def try_close_writer(writer: asyncio.StreamWriter, client_addr: tuple[str, int]):
    try:
        writer.close()
        await writer.wait_closed()
    except BrokenPipeError as e:
        log.error(
            "Error occured when closing writer for %s: %s", client_addr, e, exc_info=e
        )


# Handle for alternative implementation of TCPClient
# class ClientHandler:
#     def __init__(
#         self,
#         queue: asyncio.Queue,
#     ) -> None:
#         self.queue = queue

#     async def __call__(
#         self, reader: asyncio.StreamReader, writer: asyncio.StreamWriter
#     ):
#         addr = writer.get_extra_info("peername")
#         log.info("Client connected from %r", addr)

#         def client_cleanup(fut):
#             log.debug("Closing connection with client %s", addr)
#             try:
#                 fut.result()
#             except Exception as e:
#                 log.error(
#                     "Failed to cleanup resources on disconnecting client %s: %s",
#                     addr,
#                     e,
#                 )

#         # task = asyncio.create_task(client_task(self.queue, reader, writer))
#         # task.add_done_callback(client_cleanup)
#         await asyncio.gather(client_task(self.queue, reader, writer))


# Alternative implementation of TCPClient
# class SocketServer(Service):
#     queue: asyncio.Queue
#     host: str = hostname.value
#     port: int = port.value

#     async def start(self):
#         self.server = await asyncio.start_server(
#             ClientHandler(self.queue),
#             host=self.host,
#             port=self.port,
#             flags=socket.SOCK_STREAM | socket.AF_INET,
#         )
#         log.info("Started socket server on host %s on port %d", self.host, self.port)
#         await self.server.serve_forever()


class SocketServer(TCPServer):
    def __init__(
        self,
        host: str = hostname.value,
        port: int = port.value,
        queue: asyncio.Queue = asyncio.Queue(),
        max_connections: int = max_connections.value,
        **kwargs,
    ):
        super().__init__(address=host, port=port, **kwargs)
        self.queue = queue
        self.max_connections = max_connections
        self._connections = 0
        self._clients: dict[
            tuple[str, int], tuple[asyncio.StreamReader, asyncio.StreamWriter]
        ] = {}

    @property
    async def num_clients(self) -> int:
        return len(self._clients)

    async def handle_client(
        self, reader: asyncio.StreamReader, writer: asyncio.StreamWriter
    ) -> None:
        """Handle a client connection. In contrast to `handle_client` in `TCPClient`,
        here we handle all exceptions ourselves.

        Parameters
        ----------
        reader : asyncio.StreamReader
            a stream reader holding the underlying transport protocol
        writer : asyncio.StreamWriter
            a stream writer holding the underlying transport protocol

        Returns
        -------
        None
        """
        await asyncio.gather(self.client_task(reader, writer))
        return None

    async def client_task(
        self,
        reader: asyncio.StreamReader,
        writer: asyncio.StreamWriter,
    ):
        client_addr = writer.get_extra_info("peername")
        if is_allowed(*client_addr):
            log.debug(
                "client is connecting from %s; client is in the whitelist", client_addr
            )
            async with asyncio.Lock():
                if self._connections >= max_connections.value:
                    log.debug(
                        "Maximum number of connections reached: %d/%d",
                        self._connections,
                        max_connections.value,
                    )
                    log.debug("Client will not be connected")
                    await try_close_writer(writer, client_addr)
                    return

            await self.add_client(client_addr, reader, writer)

        else:
            log.debug(
                "client is connecting from %s; client is not allowed", client_addr
            )
            await try_close_writer(writer, client_addr)

    async def add_client(
        self,
        client_addr: tuple[str, int],
        reader: asyncio.StreamReader,
        writer: asyncio.StreamWriter,
    ):
        async with asyncio.Lock():
            self._clients[client_addr] = (reader, writer)
            log.debug("Added client %s:%d", *client_addr)

    async def remove_client(self, client_addr: tuple[str, int]):
        async with asyncio.Lock():
            removed_client = self._clients.pop(client_addr, None)
            if removed_client is not None:
                log.debug("Removed client %s:%d", *client_addr)
            else:
                log.debug("Client %s:%d not found; nothing to remove", *client_addr)

    async def watch_clients(self):
        """
        Watch for clients and remove unconnected ones.
        """
        while True:
            dead_clients = set()

            for addr, streams in self._clients.items():
                async with asyncio.Lock():
                    reader, writer = streams
                    if writer.is_closing():
                        dead_clients.add(addr)
                        log.debug("Found inactive client: %s", addr)

            for addr in dead_clients:
                await self.remove_client(addr)
                log.debug("Client removed due to inactivity: %s", addr)

            # Sleep here or the event loop will stuck
            await asyncio.sleep(0.0)

    async def send_data(self):
        """
        Send data to clients over socket connection.
        """
        while True:
            try:
                data = await self.queue.get()
                self.queue.task_done()
                if data:
                    msg = await form_message(data)

                    for client_addr, streams in self._clients.items():
                        reader, writer = streams
                        if msg:
                            log.debug("Sending data to %s", client_addr)
                            log.debug("Data: %s", msg)
                            try:
                                writer.write(msg)
                                await writer.drain()
                            except Exception as e:
                                log.debug("Client %s error: %s", client_addr, e)

                                await try_close_writer(writer, client_addr)
                                await self.remove_client(client_addr)

            except BaseException as e:
                log.error("Error sending observation plan: %s", e)
            finally:
                await asyncio.sleep(0.0)

    async def start(self) -> None:
        await asyncio.gather(super().start(), self.watch_clients(), self.send_data())
