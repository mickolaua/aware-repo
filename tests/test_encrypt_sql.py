from aware.sql.models import Subscriber, create_session, Base
from sqlalchemy.engine import create_engine, url
from sqlalchemy.orm import sessionmaker
from sqlalchemy import text
import sqlcipher3


def test():
    DB_NAME = ":memory:"
    KEY = "pass"
    engine, session = create_session(
        database=DB_NAME,
        driver="sqlite",
        passwd=KEY
    )
    Base.metadata.create_all(engine)

    with session as s:
        subs = s.query(Subscriber).all()

if __name__ == "__main__":
    test()
