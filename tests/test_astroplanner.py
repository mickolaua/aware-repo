from __future__ import annotations

from astropy.coordinates import SkyCoord, AltAz
from astropy.time import Time
from astroplan import FixedTarget
from aware.site import Telescopes


def test():
    mondy = Telescopes["mondy_azt33ik"]
    target1 = FixedTarget(
        name="GBM Trigger 683571622", coord=SkyCoord(167.3, 1.01, unit=["deg", "deg"])
    )
    target2 = FixedTarget(
        name="GBM Trigger 683644301", coord=SkyCoord(220.46, 17.23, unit=["deg", "deg"])
    )
    targets = [target1, target2]

    time1 = Time("2022-12-10T20:00:00", format="isot")
    # altaz11 = mondy.altaz(time1, target1)
    # altaz12 = mondy.altaz(time1, target2)
    # print(altaz11.secz, altaz12.secz)
    # ordered_targets = (
    #     [target1, target2] if altaz12.secz > altaz11.secz else [target2, target1]
    # )
    ordered_targets = [target2, target1]
    sorted_targets1 = mondy.observation_order(targets, time1, method="nearest")
    assert sorted_targets1 == ordered_targets, "Targets sorted incorrectly"

    time2 = Time("2022-12-10T23:59:00", format="isot")
    # altaz21 = mondy.altaz(time2, target1)
    # altaz22 = mondy.altaz(time2, target2)
    # print(altaz21.secz, altaz22.secz)
    # ordered_targets = (
    #     [target1, target2] if altaz22.secz > altaz21.secz else [target2, target1]
    # )
    ordered_targets = [target1, target2]
    # print(ordered_targets)
    sorted_targets2 = mondy.observation_order(targets, time2, method="nearest")
    # print(sorted_targets2)
    assert sorted_targets2 == ordered_targets, "Targets sorted incorrectly"


if __name__ == "__main__":
    test()
