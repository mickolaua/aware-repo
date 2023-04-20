"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
observation.py (c) 2023
Desc: a single target observation
Created:  2023-04-18
Modified: !date!
"""
from __future__ import annotations

from dataclasses import dataclass

from astroplan.target import FixedTarget
from astropy import units as u

from ..glade import GladeGalaxy

__all__ = ["Observation"]


@dataclass
class Observation:
    site_id: str = ""
    comment: str = ""
    exposure: u.Unit = 1 * u.s
    exp_count: int = 1
    filter: str = ""
    target: FixedTarget | GladeGalaxy | None = None
    sent_to: str = ""
    