"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
program.py (c) 2023
Desc: observational program for observers
Created:  2023-04-18
Modified: 2024-04-02
"""
from __future__ import annotations
import os
from typing import Any, Callable, Protocol, Sequence, Generic, TypeVar

from astropy.time import Time
from astropy import units as u
from astroplan.target import FixedTarget
from astropy.coordinates import SkyCoord

from aware.field import Field

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
        "0.00 "
        + f"{obs.exp_count:d}x{obs.exposure.to_value(u.s):.1f}"
        + (f"*{obs.filter}" if (obs.filter and obs.filter.lower() != "clear") else "")
    )

    return entry


def txt_list_entry(obs: Observation) -> str:
    ra, dec = coord2str(obs.target.coord)
    f = obs.filter if obs.filter else "clear"
    res = f"{ra} {dec} {obs.exposure.to_value(u.s):.1f} {obs.exp_count:d} {f}"
    return res


ObjectListFactory = Callable[[Sequence[Observation], bool], str]


def create_objlist(observations: Sequence[Observation], header: bool) -> str:
    return "\n".join(objlist_entry(obs) for obs in observations)


def create_txt_list(observations: Sequence[Observation], header: bool) -> str:
    hdr = "ra dec exp nframes filter\n" if header else ""
    return hdr + "\n".join([txt_list_entry(obs) for obs in observations])


ObjectListFormats: dict[str, ObjectListFactory] = {
    "objlist": create_objlist,
    "txt": create_txt_list,
}


def file_extension_to_list_format(ext: str) -> str:
    if ext in {"plan", "list", "objlist"}:
        return "objlist"
    elif ext in {"txt", "ascii"}:
        return "txt"
    else:
        raise ValueError(f"`{ext}` not recognized as a valid object plan extension")


def create_observation_program(
    site: Site, targets: Sequence[GladeGalaxy | FixedTarget]
) -> str:
    observations = [
        Observation(
            site.name,
            exposure=site.default_exposure,
            exp_count=3,
            filter=site.default_filter,
            target=tgt,
        )
        for tgt in targets
    ]
    fmt = file_extension_to_list_format(site.default_target_list_fmt)
    formatter = ObjectListFormats[fmt]

    return formatter(observations, header=site.list_formatter_kws.get("header", True))


__T = TypeVar("__T", FixedTarget, GladeGalaxy, Field)


class ProgramLoader(Generic[__T]):
    def __init__(self, target_type: __T = FixedTarget):
        self.target_type = target_type

    def __call__(self, program: str, **kwargs) -> list[__T]:
        raise NotImplementedError("__call__ must be implemented in subclass")
    

class TxtProgramLoader(ProgramLoader):
    def __call__(self, program: str, **kwargs) -> list[__T]:
        targets = []
        lines = program.splitlines()
        # Read line by line, ignoring header (first line)
        for line in lines[1:]:
            ra, dec, exptime, count, band = line.split()
            coord = SkyCoord(ra, dec, unit=["hr", "deg"])
            tgt = self.target_type(coord, **kwargs)
            targets.append(tgt)
        return targets
    

class ObjlistProgramLoader(ProgramLoader):
    def __call__(self, program: str, **kwargs) -> list[__T]:
        targets = []
        lines = program.splitlines()
        for line in lines:
            name, other_props = line.split(" = ")
            ttype, ra_str, dec_str, mag, exp_and_filter = other_props.split()
            ra = f"{ra_str[0:2]}:{ra_str[3:4]}:{ra_str[5:]}"
            dec = f"{dec_str[0:2]}:{dec_str[3:4]}:{dec_str[5:]}"
            coord = SkyCoord(ra, dec, unit=["hr", "deg"])
            tgt = self.target_type(coord, **kwargs)
            targets.append(tgt)
        return targets


def get_object_loader(ext: str) -> type[ProgramLoader]:
    fmt = file_extension_to_list_format(ext)
    if fmt == "txt":
        return TxtProgramLoader
    elif fmt == "objlist":
        return ObjlistProgramLoader
    else:
        raise NotImplementedError(f"Object loader not implemented for extension {ext}")


def read_observation_program(filename: str, object_type: __T = FixedTarget, **kwargs):
    with open(filename, 'r') as f:
        raw_program = f.read()

    basename, ext = os.path.splitext(filename)
    ext = ext[1:]
    loader = get_object_loader(ext)
    initialized_loader = loader(target_type=object_type)
    objects = initialized_loader(raw_program, **kwargs)

    return objects
    
    
