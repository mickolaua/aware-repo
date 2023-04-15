from dataclasses import dataclass, field
from datetime import datetime
from .alert.target_info import TargetInfo


@dataclass(frozen=True)
class TelegramAlertMessage:
    id: str
    senter: str
    body: str
    created: datetime = datetime.now()


@dataclass(frozen=True)
class TelegramDataPackage:
    id: str
    target_info: TargetInfo
    json_filename: str
    plot_fnames: list[str] = field(default_factory=list)
    created: datetime = datetime.now()