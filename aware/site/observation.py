from __future__ import annotations

from dataclasses import dataclass, field

from astroplan.target import FixedTarget

from ..glade import GladeGalaxy


__all__ = ["ObservationProgram"]


@dataclass
class ObservationProgram:
    site_id: str = ""
    format: str = "txt"
    comment: str = ""
    exposures: list[float] = field(default_factory=list)
    filters: list[str] = field(default_factory=list)
    targets: list[FixedTarget | GladeGalaxy] = field(default_factory=list)
