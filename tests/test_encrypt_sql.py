from aware.sql.models import Subscriber, dbconnect


def test():
    # Database keyword is not working here, due to monkey-patch in the conftest
    with dbconnect(passwd="pass") as s:
        subs = s.query(Subscriber).all()

if __name__ == "__main__":
    test()
