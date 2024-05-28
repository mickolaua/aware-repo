from aware.sql.models import Subscriber, dbconnect


def test():
    DB_NAME = ":memory:"

    subs = []
    with dbconnect(
        database=DB_NAME,
        driver="sqlite",
    ) as s:
        for i in range(0, 10):
            sub = Subscriber()
            sub.alert_type = "gcn.classic.voevent.LVC_INITIAL"
            sub.chat_id = i
            sub.telescopes = "telescope1,telescope2"
            sub.content_type = "Alerts,Schedules"
            s.add(sub)

        result = s.query(Subscriber).all()
        for res in result:
            subs.append(
                {
                    "alert_type": res.alert_type,
                    "chat_id": res.chat_id,
                    "telescopes": res.telescopes,
                    "content_type": res.content_type,
                }
            )

    assert all(
        [isinstance(sub["chat_id"], (type(None), int)) for sub in subs]
    ), "Chat_id must be integer"


if __name__ == "__main__":
    test()
