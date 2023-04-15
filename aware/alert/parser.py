from __future__ import annotations

from io import BytesIO, StringIO
from typing import Any, Protocol

from .target_info import TargetInfo


__all__ = ["AlertParser"]


# Registry for alert parser by topic
AlertParsers: dict[str, AlertParser] = {}


class AlertParser:
    topic: str = ""
    instrument: str = ""

    def __init_subclass__(cls) -> None:
        global AlertParsers
        AlertParsers[cls.topic] = cls
        return super().__init_subclass__()

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        """Parses the alert and retrieves the TargetInfo, i.e. localization
        parameters, which target it is, etc.

        Parameters
        ----------
        alert_msg : str | BytesIO | StringIO
            an alert message in form of XML-file (VOEvent)

        Returns
        -------
        TargetInfo | None
            the target information, including localization parameters,
            event name, etc.

        Raises
        ------
        NotImplementedError
            should be implemented in subclass
        """
        raise NotImplementedError

from .plugins import *
