from sqlalchemy import BLOB, Integer, REAL, String, Column, DateTime, CHAR
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.engine import create_engine, URL
from sqlalchemy.orm import sessionmaker
from .config import cfg

__all__ = ["Base", "Alert", "create_session"]

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


def create_session(database=cfg["sql"]["alert_db"]):
    url = URL.create("sqlite", database=database)
    engine = create_engine(url)

    Session = sessionmaker(bind = engine)
    session = Session()

    return engine, session
