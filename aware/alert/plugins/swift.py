from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO, StringIO
from typing import Any, Optional, Protocol

import pendulum
import requests
from astropy import units as u
from astropy.io import fits
from astropy.time import Time

from ...localization import CircularSkyMap
from ...logger import log
from ...voevent import VOEvent
from ..parser import AlertParser
from ..target_info import TargetInfo


__all__ = [
    "SWIFTXRTAlertParser",
    "SWIFTBATAlertParser",
    "SWIFT_ACTUAL_POINTDIR_AlertParser",
    "SWIFTUVOTAlertParser",
]


def describe_swift_meta(meta: dict[str, Any]) -> str:
    if meta.get("StarTrack_Lost_Lock", False):
        descr = "Star tracker lost lock\n"
    else:
        descr = (
            (f"Event also known as {meta['Name']}\n" if meta.get("Name", False) else "")
            + ("Event located near bright star\n" if meta["Near_Bright_Star"] else "")
            + (
                "Error circle located in a galaxy\n"
                if meta.get("Err_Circle_in_Galaxy", False)
                else ""
            )
            + (
                "Localization region contains a galaxy\n"
                if meta.get("Galaxy_in_Err_Circle", False)
                else ""
            )
        )

    return descr.strip("\n")


def parse_swift_alert(msg: str, topic: str, instr: str):
    # VOEvent document tree
    voevent = VOEvent.from_string(msg)
    root = voevent.root

    # Metadata
    meta = {}
    meta["GRB_Identified"] = voevent.get_parameter_value("GRB_Identified", bool)
    meta["StarTrack_Lost_Lock"] = voevent.get_parameter_value(
        "StarTrack_Lost_Lock", bool
    )
    meta["Near_Bright_Star"] = voevent.get_parameter_value("Near_Bright_Star", bool)
    meta["Err_Circle_in_Galaxy"] = voevent.get_parameter_value(
        "Err_Circle_in_Galaxy", bool
    )
    meta["Galaxy_in_Err_Circle"] = voevent.get_parameter_value(
        "Galaxy_in_Err_Circle", bool
    )

    try:
        name = root.Why.Inference.Name
    except AttributeError:
        name = ""

    meta["Name"] = name

    # We are only interested in real GRBs
    def_not_grb = voevent.get_parameter_value("Def_NOT_a_GRB", bool)
    rejected = def_not_grb or meta["StarTrack_Lost_Lock"] or not meta["GRB_Identified"]

    # Astronomical coordinates and time
    obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
    astro_coords = obs_loc.AstroCoords

    # Position of the center of the localization area
    pos2d = astro_coords.Position2D
    ra = float(pos2d.Value2.C1)
    dec = float(pos2d.Value2.C2)

    # Position uncertainty (radius)
    error_radius = float(pos2d.Error2Radius)

    # Trigger time
    isot = Time(astro_coords.Time.TimeInstant.ISOTime.text, format="isot").datetime

    # Trigger ID
    target = voevent.get_parameter_value("TrigID")

    # Trigger significance
    importance = voevent.get_parameter_value("importance", float)

    # Localization skymap
    sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)

    # Alert text description
    description = "\n".join([sky_map.describe(), describe_swift_meta(meta)])

    # Pack target info
    info = TargetInfo(
        sky_map,
        packet_type=topic,
        event=target,
        origin=instr,
        trigger_date=isot,
        importance=importance,
        meta=meta,
        description=description,
        rejected=rejected
    )

    return info


class SWIFTBATAlertParser(AlertParser):
    topic = "gcn.classic.voevent.SWIFT_BAT_GRB_POS_ACK"
    instrument = "Swift BAT"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_swift_alert(
            alert_msg, SWIFTBATAlertParser.topic, SWIFTBATAlertParser.instrument
        )


class SWIFTXRTAlertParser(AlertParser):
    topic = "gcn.classic.voevent.SWIFT_XRT_POSITION"
    instrument = "Swift XRT"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_swift_alert(
            alert_msg, SWIFTXRTAlertParser.topic, SWIFTXRTAlertParser.instrument
        )


class SWIFTUVOTAlertParser(AlertParser):
    topic = "gcn.classic.voevent.SWIFT_UVOT_POS"
    instrument = "Swift UVOT"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_swift_alert(
            alert_msg, SWIFTUVOTAlertParser.topic, SWIFTUVOTAlertParser.instrument
        )


class SWIFT_ACTUAL_POINTDIR_AlertParser(AlertParser):
    topic = "gcn.classic.voevent.SWIFT_ACTUAL_POINTDIR"
    instrument = "Swift BAT"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_swift_alert(
            alert_msg,
            SWIFT_ACTUAL_POINTDIR_AlertParser.topic,
            SWIFT_ACTUAL_POINTDIR_AlertParser.instrument,
        )
