"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_avro.py (c) 2024
Desc: testing parsing of avro schema
Created:  2024-08-10
Modified: !date!
"""
import orjson
from aware.alert.plugins.ep import WXTAlertParser
from astropy.time import Time
from pytest import approx


def test():
    ALERT = {
        "$schema": "https://gcn.nasa.gov/schema/v4.1.0/gcn/notices/einstein_probe/wxt/"
                   "alert.schema.json",
        "instrument": "WXT",
        "trigger_time": "2024-03-01T21:46:05.13Z",
        "id": [
            "01708973486"
        ],
        "ra": 120,
        "dec": 40,
        "ra_dec_error": 0.02,
        "image_energy_range": [
            0.5,
            4
        ],
        "net_count_rate": 1,
        "image_snr": 1,
        "additional_info": "The net count rate is derived from an accumulated image "
                           "(up to 20 min) in 0.5-4 keV, assuming a constant flux. "
                           "However, it can be significantly lower than the actual "
                           "count rate of a burst with a duration much shorter than 20 "
                           "min."
    }

    ALERT_BYTES = orjson.dumps(ALERT)
    target_info = WXTAlertParser.parse_alert(ALERT_BYTES)
    assert target_info.event == ALERT["id"][0]
    assert target_info.trigger_date == Time(ALERT["trigger_time"], format="isot")
    assert target_info.localization.center().ra.deg == ALERT["ra"]
    assert target_info.localization.center().dec.deg == ALERT["dec"]
    assert target_info.importance < 0.9
    assert target_info.localization.error_radius().to_value("deg") == approx(0.02)


if __name__ == "__main__":
    test()


