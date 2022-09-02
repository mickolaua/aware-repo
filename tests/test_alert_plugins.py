from __future__ import annotations

# from aware import plugin_manager
from aware.alert import plugin_manager
from aware.logger import log


def test():
    res = plugin_manager.hook.parse_alert(alert_msg="tests/event.xml")
    log.info(res)
    
    assert res.ra_center == 0.0, "Parser plugin returned incorrect R.A."


if __name__ == "__main__":
    test()

