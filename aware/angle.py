from __future__ import annotations

from astropy.coordinates import SkyCoord, Angle


def coord2str(coord: SkyCoord) -> tuple[str, str]:
    return coord.to_string("hmsdms", sep=":", precision=2, pad=True).split(" ")


def coord_to_target_name(coord: SkyCoord) -> str:
    ra, dec = coord.to_string("hmsdms", sep="", precision=2, pad=True).split(" ")
    return f"J{ra}_{dec}".replace(".", "_")
