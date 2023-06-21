from __future__ import annotations
import time

from typing import Any
import posixpath
from uuid import uuid4

import requests
from sqlalchemy import and_, exc, func, join, select

from ..logger import log
from ..sql.models import Alert, RejectedAlert, create_session
from ..data import products_dir


class RejectedAlertNotFound(exc.NoResultFound):
    def __init__(self, *arg: Any, event_id: str, **kw: Any):
        super().__init__(*arg, **kw)
        self.event_id = event_id

    def _message(self) -> str:
        return f"retraction information not found for event {self.event_id}"


def is_retracted(event: str, origin: str) -> bool:
    _, session = create_session()
    with session as s:
        stmnt = (
            select(func.count())
            .select_from(join(Alert, RejectedAlert, Alert.id == RejectedAlert.alert_id))
            .where(and_(Alert.event == event, Alert.origin == origin))
        )
        res = s.scalar(stmnt)
        if res is not None:
            return res > 0
        else:
            return False


def is_retracted_by_id(alert_id: int) -> bool:
    _, session = create_session()
    with session as s:
        stmnt = (
            select(func.count())
            .select_from(join(Alert, RejectedAlert, Alert.id == RejectedAlert.alert_id))
            .where(Alert.id == alert_id)
        )
        res = s.scalar(stmnt)
        if res is not None:
            return res > 0
        else:
            return False


def add_retracted(alert_id: int, event: str, origin: str):
    _, session = create_session()
    with session as s:
        try:
            s.add(RejectedAlert(alert_id=alert_id, event=event, origin=origin))
            s.flush()
        except exc.IntegrityError as e:
            s.rollback()

        try:
            s.commit()
        except exc.OperationalError:
            s.rollback()


def max_id() -> int | None:
    _, session = create_session()
    with session as s:
        try:
            stmt = select(func.count()).select_from(Alert)
            return s.execute(stmt).scalar_one_or_none()
        except Exception as e:
            log.error("maximal alert_id no retrieved", exc_info=e)
            return None

