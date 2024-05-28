from __future__ import annotations

from pathlib import Path
from typing import Sequence

from aware.glade import GladeGalaxy
from .alert.target_info import TargetInfo
from astroplan.target import FixedTarget
from .site import Site
import json


class JSON:
    def __init__(
        self, 
        targets: Sequence[GladeGalaxy | FixedTarget], 
        site: Site, 
        info: TargetInfo, comment: str = ""
    ) -> None:
        self._targets = targets
        self._site = site
        self._info = info
        self._comment = comment

        self.json = {}
        self.json["site"] = {
            "name": self._site.name, 
            "lon": self._site.location.lon.value, 
            "lat": self._site.location.lat.value,
        }
        center = self._info.localization.center()
        radius = self._info.localization.error_radius().to_value("deg")
        self.json["event"] = {
            "name": self._info.event,
            "ra": center.ra.deg,
            "dec": center.dec.deg,
            "error_radius": radius,
        }
        self.json["targets"] = [
            {   
                "name": str(t.name),
                "ra": t.coord.ra.value, 
                "dec": t.coord.dec.value,
                "Bmag": getattr(t, "Bmag", None),
                "z": getattr(t, "z", None),
            }
            for t in self._targets
        ]
        self.json["comment"] = self._comment
        self._structures = [self.json]

    def add_jsons(self, *jsons: JSON) -> list[JSON]:
        self._structures += [j.json for j in jsons]

    def to_string(self) -> str:
        return json.dumps(self._structures, indent=2)
