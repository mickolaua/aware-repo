"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
ep.py (c) 2024
Desc: Alert parser for Enstein Probe (EP)
Created:  2024-05-17
Modified: !date!
"""

from io import BytesIO, StringIO
import orjson
from astropy.time import Time
from aware.logger import log
from aware.alert.parser import AlertParser
from aware.alert.target_info import TargetInfo
from aware.localization.main import CircularSkyMap


class WXTAlertParser(AlertParser):
    topic = "gcn.notices.einstein_probe.wxt.alert"
    instrument = "EP-WXT"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args, **kwargs
    ) -> TargetInfo | None:
        try:
            msg_content = orjson.loads(alert_msg)
        except orjson.JSONDecodeError as e:
            log.error("Could not parse alert due to incorrect JSON encoding: %s", e)
            info = None
        else:
            ra = msg_content["ra"]
            dec = msg_content["dec"]
            radius = msg_content["ra_dec_error"]
            trigger_date = Time(
                msg_content["trigger_time"], format="isot"
            ).to_datetime()
            importance = 0.9
            event = f"J{ra:.3f}{dec:+.3f}"
            localization = CircularSkyMap(ra_center=ra, dec_center=dec, radius=radius)
            info = TargetInfo(
                localization=localization,
                packet_type=WXTAlertParser.topic,
                event=event,
                origin=WXTAlertParser.instrument,
                trigger_date=trigger_date,
                importance=importance,
            )
        return info
