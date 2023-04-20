from __future__ import annotations
from dataclasses import dataclass

from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord


@dataclass
class Box:
    left_lower: SkyCoord
    right_lower: SkyCoord
    left_upper: SkyCoord
    right_upper: SkyCoord


class Field(FixedTarget):
    def __init__(
        self,
        center_coord: SkyCoord,
        width: Angle = Angle(1 * u.deg),
        height: Angle = Angle(1 * u.deg),
        name: str | None = None,
        **kwargs,
    ):
        super().__init__(center_coord, name, **kwargs)
        self.width = width
        self.height = height

    def box(self) -> Box:
        """
        Return coordinates of the field corners.
        """
        # Left upper corner
        ra_left = self.coord.ra - self.width*0.5
        dec_top = self.coord.dec + self.height*0.5
        left_upper = SkyCoord(ra_left, dec_top)

        # Right upper corner
        ra_right = self.coord.ra + self.width*0.5
        right_upper = SkyCoord(ra_right, dec_top)

        # Left lower corner
        dec_bottom = self.coord.dec - self.width*0.5
        left_lower = SkyCoord(ra_left, dec_bottom)

        # Right lower corner
        right_lower = SkyCoord(ra_right, dec_bottom)

        return Box(left_lower, right_lower, left_upper, right_upper)


    def __str__(self) -> str:
        return super().__str__()
