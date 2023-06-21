from __future__ import annotations
import logging

import os.path
import pytest
from sqlalchemy.orm import sessionmaker, object_session
from sqlalchemy.engine import create_engine

import aware.sql.models
import aware.alert.crossmatch


@pytest.fixture
def test_dir():
    DIR = os.path.dirname(__file__)
    return DIR


BIND = SESSION = None


def create_session(**kwargs):
    global BIND, SESSION
    if not BIND or not SESSION:
        BIND = create_engine(url="sqlite+pysqlite:///:memory:")
        aware.sql.models.Base.metadata.create_all(BIND)
        SessionMaker = sessionmaker(BIND)
        SESSION = SessionMaker()

    def f(**kwargs):
        return BIND, SESSION

    return f


@pytest.fixture(autouse=True)
def mocksession1(monkeypatch):
    monkeypatch.setattr(aware.sql.models, "create_session", create_session())
    

@pytest.fixture(autouse=True)
def mocksession2(monkeypatch):
    monkeypatch.setattr(
        aware.alert.crossmatch, "create_session", create_session()
    )


@pytest.fixture(autouse=True)
def no_logging():
    logger = logging.getLogger("aware")
    logger.disabled = True
