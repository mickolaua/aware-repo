"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
program.py (c) 2023
Desc: observational program for observers
Created:  2023-04-18
Modified: !date!
"""
from __future__ import annotations
from typing import Callable, Sequence

from astropy.time import Time
from astropy import units as u
from astroplan.target import FixedTarget

from .observation import Observation
from ..glade import GladeGalaxy
from ..angle import coord2str
from ..site import Site


def objlist_entry(obs: Observation) -> str:
    """Create an entry of the .objlist observational program for CHAOS/FORTE software.

    Parameters
    ----------
    obs : Observation
        an observation of the certain target (e.g. Glade+ galaxy or sky field)

    Returns
    -------
    entry: str
        string representation of the .objlist entry
    """
    ra, dec = coord2str(obs.target.coord)
    ra = ra.replace(":", "")
    dec = dec.replace(":", "")
    
    entry = (
        f"{obs.target.name} = F {ra} {dec} "
        +f"{obs.exp_count}x{obs.exposure.to_value(u.s):.1f}*"
        +f"{obs.filter} !10"
    )

    return entry


def txt_list_entry(obs: Observation) -> str:
    ra, dec = coord2str(obs.target.coord)
    res = f"{ra} {dec} {obs.exposure.to_value(u.s)}x{obs.exp_count} 10"
    return res


ObjectListFactory = Callable[[Sequence[Observation]], str]


def create_objlist(observations: Sequence[Observation]) -> str:
    return "\n".join(objlist_entry(obs) for obs in observations)


def create_txt_list(observations: Sequence[Observation]) -> str:
    hdr = "ra dec exp filter priority\n"
    return hdr + "\n".join([txt_list_entry(obs) for obs in observations])


ObjectListFormats: dict[str, ObjectListFactory] = {
    "objlist": create_objlist,
    "txt": create_txt_list
}


def create_observation_program(
    site: Site, targets: Sequence[GladeGalaxy | FixedTarget]
) -> str:
    observations = [
        Observation(
            site.name, 
            exposure=site.default_exposure, 
            exp_count=3, 
            filter=site.default_filter, 
            target=tgt
        )
        for tgt in targets
    ]

    return ObjectListFormats[site.default_target_list_fmt](observations)