from aware.sql.models import Subscriber, create_session, Base


def test():
    DB_NAME = ":memory:"
    engine, session = create_session(
        database=DB_NAME,
        driver="sqlite",
    )
    Base.metadata.create_all(engine)

    with session as s:
        subs = []
        for i in range(0, 10):
            sub = Subscriber()
            sub.alert_type = "gcn.classic.voevent.LVC_INITIAL"
            sub.chat_id = i
            sub.telescopes = "telescope1,telescope2"
            sub.content_type = "Alerts,Schedules"
            subs.append(sub)
        s.add_all(subs)
        s.commit()

        subs = s.query(Subscriber).all()
        assert all(
            [isinstance(sub.chat_id, (type(None), int)) for sub in subs]
        ), "Chat_id must be integer"


if __name__ == "__main__":
    test()
