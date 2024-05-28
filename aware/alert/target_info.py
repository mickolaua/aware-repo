from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Optional

from ..localization import Localization


__all__ = ["TargetInfo"]


@dataclass
class TargetInfo:
    localization: Localization | None
    packet_type: str | None = ""
    event: str | None = ""
    origin: str | None = ""
    trigger_date: Optional[datetime] = None
    importance: Optional[float] = None
    description: Optional[str] = None
    rejected: bool | None = False
    meta: Optional[dict[Any, Any]] = field(default_factory=dict)

    def describe(self) -> str:
        return self.localization.describe()