"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_alert_crossmatch.py (c) 2023
Desc: test the crossmatching of the alert pair
Created:  2023-05-07
Modified: !date!
"""
from __future__ import annotations
import asyncio
import logging

from aware.alert import AlertParsers
from aware import alert
from aware.sql.models import Base, RejectedAlert
import os.path
from aware.sql.models import Alert, create_session
from aware.alert.crossmatch import crossmatch_alerts, replace_with_matched
from astropy.time import Time
import pytest
from aware.consumer.main import dump_alert_to_db
from unittest import mock
from sqlalchemy.orm import sessionmaker
from sqlalchemy.engine import create_engine
from aware.sql.models import Base
import aware.sql.models
import aware.alert.crossmatch 

logger = logging.getLogger("aware")
logger.disabled = True


@pytest.mark.asyncio
async def test(test_dir):
    parser = AlertParsers['gcn.classic.voevent.SWIFT_XRT_POSITION']
    
    DIR = test_dir
    with open(os.path.join(DIR, "alert_messages/xrt_alert_1165948_acc.xml"), "rb") as f:
        alert_msg_acc = f.read()

    with open(os.path.join(DIR, "alert_messages/xrt_alert_1165948.xml"), "rb") as f:
        alert_msg = f.read()

    info_acc = parser.parse_alert(alert_msg_acc)
    info_acc.rejected = False

    info = parser.parse_alert(alert_msg)
    info.rejected = False

    await dump_alert_to_db(alert_msg_acc, info_acc)
    # await dump_alert_to_db(alert_msg, info)

    m = await crossmatch_alerts(info)
    assert len(m) == 1, "There is only one matched event!"

    # if info.rejected:
    #     for al in m:
    #         alert.util.add_retracted(al.id, al.event, al.origin)

    # with session as s:
    #     rej_alerts = s.query(RejectedAlert).all()

    #     for al in rej_alerts:
    #         print(al.alert.id)

    # await replace_with_matched(m, alert_msg_acc, info_acc)

    # _, session = create_session()
    # with session as s:
    #     alert_acc = Alert()
    #     alert_acc.alert_message = alert_msg_acc
    #     alert_acc.dec_center = info_acc.localization.center().dec.deg
    #     alert_acc.error_radius = info_acc.localization.error_radius().to_value("deg")
    #     alert_acc.event = info_acc.event
    #     alert_acc.importance = info_acc.importance
    #     alert_acc.localization = bytes(str(info_acc.localization), encoding="utf-8")
    #     alert_acc.origin = info_acc.origin
    #     alert_acc.ra_center = info_acc.localization.center().ra.deg
    #     alert_acc.trigger_date = info_acc.trigger_date

    #     alert = Alert()
    #     alert.alert_message = alert_msg_acc
    #     alert.dec_center = info.localization.center().dec.deg
    #     alert.error_radius = info.localization.error_radius().to_value("deg")
    #     alert.event = info.event
    #     alert.importance = info.importance
    #     alert.localization =  bytes(str(info.localization), encoding="utf-8")
    #     alert.origin = info.origin
    #     alert.ra_center = info.localization.center().ra.deg
    #     alert.trigger_date = info.trigger_date

    #     s.add_all([alert, alert_acc])
    #     s.commit()


if __name__ == "__main__":
    asyncio.run(test())
