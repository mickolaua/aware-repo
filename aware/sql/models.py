from __future__ import annotations

from typing import Any

import pandas as pd
from sqlalchemy import BLOB, CHAR, REAL, Column, DateTime, Integer, String
from sqlalchemy.engine import URL, Engine, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from ..config import CfgOption

__all__ = ["Base", "Alert", "create_session"]


alert_db = CfgOption("alert_db", "alert.db", str)
Base = declarative_base()


class Alert(Base):
    __tablename__ = "alert"

    alert_id = Column(Integer, autoincrement=True, primary_key=True)
    alert_message = Column(BLOB)
    ra_center = Column(REAL)
    dec_center = Column(REAL)
    error_radius1 = Column(REAL)
    error_radius2 = Column(REAL)
    localization = Column(BLOB)
    trigger_date = Column(DateTime)
    event = Column(String(256))
    origin = Column(String(256))
    importance = Column(REAL)
    md5 = Column(CHAR(32), nullable=False, unique=True)


def create_session(database: str = alert_db.get_value()) -> tuple[Engine, Any]:
    url = URL.create("sqlite", database=database)
    engine = create_engine(url)

    Session = sessionmaker(bind=engine)
    session = Session()

    return engine, session


def select_from_table(
    engine, session, tablename: str, column: str, value: str
) -> pd.DataFrame | None:
    df = pd.read_sql_table(tablename, engine)
    try:
        result = df[df[column] == value]
    except (LookupError, ArithmeticError, TypeError) as e:
        result = None

    return result
