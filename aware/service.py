import asyncio
import uuid

from aiomisc import Service
from gcn_kafka.oidc import set_oauth_cb

from .consumer import ConsumeLoop, prepare_consumer
from .credentials import get_credentials
from .socket import SocketServer
from .telegram.bot import TelegramBot


class ServiceFactory:
    COMMON_CONSUMER_OPTIONS = {
        "group.id": uuid.uuid4().hex,
        # Infinitely read messages in development mode since the beginning of
        # the partition
        "enable.auto.commit": False,
        # Set maximum idle time out between reading topic records to 1 day
        # (in msec) to prevent application maximum poll interval exceeded
        "max.poll.interval.ms": 86_400_000,
        "sasl.oauthbearer.config": "oauth_cb",
        "oauth_cb": set_oauth_cb,
    }

    def __init__(self, queue: asyncio.Queue):
        self.queue = queue

    def create_test_consumer(self, queue: asyncio.Queue = None) -> Service:
        cred = get_credentials()
        test_consumer_config = self.COMMON_CONSUMER_OPTIONS.copy()
        test_consumer_config.update(
            {
                "auto.offset.reset": "earliest",
                "enable.auto.commit": False,
                "debug": "generic,security,broker",
            }
        )
        consumer = prepare_consumer(test_consumer_config, cred)
        return ConsumeLoop(consumer, queue or self.queue)

    def create_prod_consumer(self, queue: asyncio.Queue = None):
        cred = get_credentials()
        prod_consumer_config = self.COMMON_CONSUMER_OPTIONS.copy()
        prod_consumer_config.update(
            {
                "auto.offset.reset": "latest",
            }
        )
        consumer = prepare_consumer(prod_consumer_config, cred)
        return ConsumeLoop(consumer, queue or self.queue)

    def create_telegram_bot(self, queue: asyncio.Queue = None) -> Service:
        bot = TelegramBot(queue=queue or self.queue)
        return bot

    def create_socket_server(self, queue: asyncio.Queue = None) -> Service:
        server = SocketServer(queue=queue or self.queue)
        return server
