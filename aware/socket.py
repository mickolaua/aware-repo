import asyncio
import os

import aiofiles
from aiomisc.service import TCPServer

from aware.config import CfgOption
from aware.data import AlertMessage, DataPackage
from aware.logger import log

hostname = CfgOption("hostname", "127.0.0.1", str)
port = CfgOption("port", 55555, int)
max_connections = CfgOption("max_connections", 5, int)
send_alert_message = CfgOption("send_alert_message", False, bool)
send_cancelled_alerts = CfgOption("send_cancelled_alerts", False, bool)


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


async def client_task(
    queue: asyncio.Queue, reader: asyncio.StreamReader, writer: asyncio.StreamWriter
):
    client_addr = writer.get_extra_info("peername")
    while True:
        try:
            data = await queue.get()

            if data:
                msg = await form_message(data)
                if msg:
                    log.debug("Sending observation plan to %s", client_addr)
                    log.debug("Plan: %s", msg)
                    writer.write(msg)
                    await writer.drain()

            queue.task_done()

        except BaseException as e:
            log.error("Error sending observation plan: %s", e)
            break

    try:
        writer.close()
        await writer.wait_closed()
    except BrokenPipeError as e:
        log.error("Error occured at closing writer for %s: %s", client_addr, e)


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
        **kwargs,
    ):
        super().__init__(address=host, port=port, **kwargs)
        self.queue = queue

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
        await asyncio.gather(client_task(self.queue, reader, writer))
        return None
