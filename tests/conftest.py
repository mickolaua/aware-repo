from __future__ import annotations

import logging
import os.path

import pytest
from sqlalchemy.engine import create_engine
from sqlalchemy.orm import sessionmaker

import aware.alert.crossmatch
import aware.sql.models


@pytest.fixture
def test_dir():
    DIR = os.path.dirname(__file__)
    return DIR


BIND = SESSION = None


def create_session(**kwargs):
    global BIND, SESSION
    if not BIND or not SESSION:
        BIND = create_engine(url=kwargs["database"])
        aware.sql.models.Base.metadata.create_all(BIND)
        SessionMaker = sessionmaker(BIND)
        SESSION = SessionMaker()

    def f(**kwargs):
        return BIND, SESSION

    return f


@pytest.fixture(autouse=True)
def mocksession1(monkeypatch, tmp_path):
    DB_DIR = tmp_path
    DB_DIR.mkdir(exist_ok=True)
    DB_NAME = "test.db"
    DB_URL = f"sqlite+pysqlite:///{DB_DIR}/{DB_NAME}"
    monkeypatch.setattr(
        aware.sql.models,
        "create_session",
        create_session(database=DB_URL, driver="sqlite"),
    )


@pytest.fixture(autouse=True)
def no_logging():
    logger = logging.getLogger("aware")
    logger.disabled = True
