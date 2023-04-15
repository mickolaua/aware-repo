from dataclasses import dataclass
import os


@dataclass(frozen=True)
class Credentials:
    id: str
    secret: str


def get_credentials() -> Credentials:
    # Credentials
    try:
        CLIENT_ID = os.environ["GCN_KAFKA_CLIENT_ID"]
    except KeyError as e:
        raise RuntimeError("client id for GCN Kafka broker is not set")

    try:
        CLIENT_SECRET = os.environ["GCN_KAFKA_CLIENT_SECRET"]
    except KeyError as e:
        raise RuntimeError(
            "client secret for GCN Kafka broker is not specified"
        )

    return Credentials(id=CLIENT_ID, secret=CLIENT_SECRET)
