from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Optional

from ..localization import Localization


__all__ = ["TargetInfo"]


@dataclass
class TargetInfo:
    localization: Localization
    packet_type: str = ""
    event: str = ""
    origin: str = ""
    trigger_date: Optional[datetime] = None
    importance: Optional[float] = None
    description: Optional[str] = None
    meta: Optional[dict[Any, Any]] = field(default_factory=dict)

    def describe(self) -> str:
        return self.localization.describe()