from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Generator, Protocol

import sqlcipher3
from sqlalchemy import (
    BLOB,
    REAL,
    VARCHAR,
    Column,
    DateTime,
    Dialect,
    ForeignKey,
    Integer,
    MetaData,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.engine import URL, Engine, create_engine
from sqlalchemy.orm import (
    DeclarativeBase,
    Session,
    relationship,
    sessionmaker,
    scoped_session,
)
from sqlalchemy.sql.operators import OperatorType
from sqlalchemy.types import TypeDecorator

from ..config import CfgOption
from aware.logger import log

__all__ = ["Base", "Alert", "Subscriber", "create_session"]


alert_db = CfgOption("alert_db", "alert.db", str, comment="Path to alert database")
driver_db = CfgOption("driver", "sqlite", str, comment="Driver to use with the database (e.g. sqlite for SQLite databases)")
user_db = CfgOption("user", "", str, comment="Database user")
password_db = CfgOption("password", "", str, comment="Database password for the user")
host_db = CfgOption("host", "", str, comment="Database hostname")
port_db = CfgOption("port", 8080, int, comment="Database port number")
query_db = CfgOption("query", "", list, comment="Database queries to be executed at the first connection")


naming_convention = {
    "ix": "ix_%(column_0_label)s",
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(column_0_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s",
}
metadata = MetaData(naming_convention=naming_convention)


class Base(DeclarativeBase):
    metadata = metadata


class TelegramID(TypeDecorator):
    """
    Telegram ID, which is stored as a string in the database, but is returned as a
    integer in Python.
    """

    impl = VARCHAR
    cache_ok = True

    def process_bind_param(self, value: Any | None, dialect: Dialect) -> Any:
        if value is not None:
            value = str(value)
        return value

    def process_result_value(self, value: Any | None, dialect: Dialect) -> Any | None:
        if value is not None:
            value = int(value)
        return value
    
    def coerce_compared_value(self, op: OperatorType | None, value: Any) -> Any:
        return VARCHAR(256)


class Alert(Base):
    __tablename__ = "alert"

    id = Column(Integer, autoincrement=True, primary_key=True)
    alert_message = Column(Text)
    ra_center = Column(REAL)
    dec_center = Column(REAL)
    error_radius = Column(REAL)
    localization = Column(BLOB)
    trigger_date = Column(DateTime)
    event = Column(String(256))
    origin = Column(String(256))
    importance = Column(REAL)


class RejectedAlert(Base):
    __tablename__ = "reject_alert"

    id = Column(Integer, autoincrement=True, primary_key=True)
    event = Column(String(256))
    origin = Column(String(256))

    # One-to-one relationship with alert table.
    alert_id = Column(Integer, ForeignKey(Alert.id))
    alert = relationship(Alert, backref=Alert.__tablename__)


class MatchedAlert(Base):
    __tablename__ = "matched_alert"
    id = Column(Integer, primary_key=True)
    alert_id = Column(Integer, ForeignKey(Alert.id), unique=True)
    event = Column(String(256))
    origin = Column(String(256))


class Subscriber(Base):
    __tablename__ = "settings"
    id = Column(Integer, autoincrement=True, primary_key=True)
    chat_id = Column(TelegramID(256), unique=True)
    content_type = Column(String(256))
    alert_type = Column(Text)
    telescopes = Column(Text)


# class FileMessage(Base):
#     __tablename__ = "file_message"
#     id = Column(Integer, autoincrement=True, primary_key=True)
#     filename = Column(Text)
#     file_id = Column(TelegramID, unique=True)


def _create_sqlite_encrypted_url(
    password: str,
    database: str,
    algorithm: str = "aes-256-cfb",
    params: str = "kdf_iter=64000",
) -> str:
    engine_name = "sqlite+pysqlcipher"
    url = f"{engine_name}://:{password}@/{database}?cipher={algorithm}&{params}'"
    return url


def _create_url(
    driver: str,
    database: str,
    user: str,
    passwd: str,
    host: str,
    port: int,
    query: dict[str, str],
) -> URL | str:
    # Regular SQLite does not provide login, but when it is encrypted, we can pass
    # encryption key as a password and decrypt it.
    if driver == "sqlite":
        if query:
            query = {}

        if user:
            user = None

        if port:
            port = None

        if passwd:
            return _create_sqlite_encrypted_url(passwd, database)

    # Otherwise just use the sqlalchemy built-in solution
    return URL.create(
        drivername=driver,
        database=database,
        username=user,
        password=passwd,
        host=host,
        port=port,
        query=query,
    )


def create_session(
    database: str = alert_db.get_value(),
    driver: str = driver_db.get_value(),
    user: str = user_db.get_value(),
    passwd: str = password_db.get_value(),
    host: str = host_db.get_value(),
    port: int = port_db.get_value(),
    query: dict[str, str] | list[str] = query_db.get_value(),
) -> tuple[Engine, Session]:
    # Convert query to dictionary
    if isinstance(query, list):
        keys = [kvpair.split("=")[0].rstrip() for kvpair in query]
        vals = [kvpair.split("=")[1].rstrip() for kvpair in query]
        query = dict(zip(keys, vals))

    url = _create_url(driver, database, user, passwd, host, port, query)
    options = dict(url=url, future=True)

    # if "sqlite" in driver and passwd:
    #     engine = create_engine(url, module=sqlcipher3, future=True)
    # else:
    #     engine = create_engine(url, future=True)
    if "sqlite" in driver and passwd:
        options.update(dict(module=sqlcipher3))

    if "sqlite" in driver:
        options.update(connect_args=dict(timeout=10))

    engine = create_engine(**options)

    session_factory = sessionmaker(bind=engine, future=True, expire_on_commit=False)
    Session = scoped_session(session_factory)

    return engine, Session


class SQLExceptionHandler(Protocol):
    """
    Protocol class for custom SQL session exception handlers. A user should implement
    a custom function that takes exactly one required argument -- exception value, and
    any number of additional arguments and keyword arguments. It is intended that
    handling of the different exceptions is performed via using isinstance function.
    Example usage below:
    ```
    def my_exception_handler(exception):
        if isinstance(exception, OperationError):
            log.exception("Operation error")
        elif isinstance(exception, AttributeError):
            log.exception("'Session' object has no attribute 'insert'")
    ```
    """

    def __call__(self, exc_value: BaseException, *args: Any, **kwds: Any) -> None:
        raise NotImplementedError("Implement exception handler in protocols")


def default_sql_handler(exc_value: BaseException):
    log.error("SQL session error: %s" % exc_value, exc_info=exc_value)



@contextmanager
def dbconnect(
    exception_handler: SQLExceptionHandler = default_sql_handler, **session_kws
) -> Generator[Session]:
    """
    Context manager for connecting to the database and committing changes.
    See ``create_session`` for more information on kwargs.
    """
    eng, session = create_session(**session_kws)
    try:
        yield session
    finally:
        try:
            session.commit()
        except Exception as e:
            exception_handler(e)
            session.rollback()
        finally:
            session.close()
