from __future__ import annotations
from dataclasses import dataclass

from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from .angle import coord_to_target_name


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
        self.name = (
            name
            if name
            else (
                f"Field {self.width.to_value(u.deg):.1f}d x "
                f"{self.height.to_value(u.deg):.1f}d at the {self.coord}"
            )
        )

    def box(self) -> Box:
        """
        Return coordinates of the field corners.
        """
        # Left upper corner
        ra_left = max(self.coord.ra.deg - self.width.to_value(u.deg) * 0.5, 0) * u.deg
        dec_top = (
            min(self.coord.dec.deg + self.height.to_value(u.deg) * 0.5, 90) * u.deg
        )
        left_upper = SkyCoord(ra_left, dec_top)

        # Right upper corner
        ra_right = (
            min(self.coord.ra.deg + self.width.to_value(u.deg) * 0.5, 360) * u.deg
        )
        right_upper = SkyCoord(ra_right, dec_top)

        # Left lower corner
        dec_bottom = (
            max(self.coord.dec.deg - self.width.to_value(u.deg) * 0.5, -90) * u.deg
        )
        left_lower = SkyCoord(ra_left, dec_bottom)

        # Right lower corner
        right_lower = SkyCoord(ra_right, dec_bottom)

        return Box(left_lower, right_lower, left_upper, right_upper)

    def __str__(self) -> str:
        return self.name
