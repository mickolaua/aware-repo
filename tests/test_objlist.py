from __future__ import annotations
from astroplan.target import FixedTarget
from astropy.coordinates import SkyCoord
from astropy import units as u
from aware.planning.observation import Observation
from aware.planning.program import create_observation_program
from aware.site import Telescopes


def test():
    azt = Telescopes["altai_santel400"]
    target = FixedTarget(
        coord=SkyCoord(150*u.deg, +35.5*u.deg), name="test_target"
    )
    program = create_observation_program(azt, [target])
    assert program == "test_target = F 100000.00 +353000.00 0.00 3x120.0"


if __name__ == "__main__":
    test()

