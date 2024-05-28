"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
integral.py (c) 2023
Desc: parsers of the alert messages by INTEGRAL spacecraft
Created:  2023-07-10
Modified: 2024-04-04
"""
from __future__ import annotations

from io import BytesIO, StringIO
from typing import Any

from astropy import units as u
from astropy.io import fits
from astropy.time import Time

from ...localization import CircularSkyMap
from ...logger import log
from ...voevent import VOEvent
from ..parser import AlertParser
from ..target_info import TargetInfo
from ...config import dev


__all__ = ["IntegralWeakParser", "IntegralSPIACSParser"]


INTEGRAL_SPI_ACS = "INTEGRAL (SPI-ACS)"
INTEGRAL_MANY = "INTEGRAL (at-least one instrument)" 


def parse_integral_voevent(
    alert_msg: str | BytesIO | StringIO, instrument: str = "", topic: str = ""
) -> TargetInfo | None:
    # VOEvent document tree
    event = VOEvent.from_string(alert_msg)
    root = VOEvent.from_string(alert_msg).root
    log.debug(
        "Received VOEvent from %s of type %s with body \n%s",
        instrument,
        topic,
        alert_msg,
    )

    # Check if alert message is test
    is_test = event.get_parameter_value("Test_Notice", type_=bool)
    if is_test and not dev.value:
        return None

    # If it GRB or not
    notgrb = event.get_parameter_value("Def_NOT_a_GRB", type_=bool)
    if notgrb:
        retraction = True
    else:
        retraction = False

    # Astronomical coordinates and time
    try:
        obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
    except AttributeError:
        obs_loc = None
    
    try:
        astro_coords = obs_loc.AstroCoords
    except AttributeError:
        astro_coords = None

    # Position of the center of the localization area
    try:
        pos2d = astro_coords.Position2D
        ra = float(pos2d.Value2.C1)
        dec = float(pos2d.Value2.C2)
    except AttributeError:
        pos2d = None
        ra = None
        dec = None

    # Position uncertainty (radius) statistical only
    try:
        error_radius = float(pos2d.Error2Radius)
    except AttributeError:
        error_radius = None

    # Trigger time
    try:
        isot = Time(astro_coords.Time.TimeInstant.ISOTime.text, format="isot").datetime
    except AttributeError:
        isot = None

    # Trigger ID
    target = event.get_parameter_value("TrigID")

    # Trigger significance
    try:
        importance = 1.0 - float(root.Why.Inference.attrib["probability"])
    except AttributeError:
        importance = 1.0

    # Pack target info, remember SPI-ACS does not provide coordinates and errors
    if not error_radius and not ra and not dec:
        sky_map = None
    else:
        sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)

    descr = ""
    if retraction:
        descr += "\nEVENT IS NOT A GRB"

    # Collect metadata
    meta = {
        "Burst_TJD": event.get_parameter_value("Burst_TJD", type_=float),
        "Burst_SOD": event.get_parameter_value("Burst_SOD", type_=float),
        "Inten_Sigma": event.get_parameter_value("Inten_Sigma", type_=float),
        "Time_Error": event.get_parameter_value("Time_Error", type_=float),
        "Time_Scale": event.get_parameter_value("Time_Scale", type_=float)
    }

    info = TargetInfo(
        sky_map,
        packet_type=topic,
        event=target,
        origin=instrument,
        trigger_date=isot,
        importance=importance,
        rejected=retraction,
        description=descr,
        meta=meta
    )

    return info


class IntegralWeakParser(AlertParser):
    topic = "gcn.classic.voevent.INTEGRAL_WEAK"
    instrument = INTEGRAL_MANY

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_integral_voevent(
            alert_msg, IntegralWeakParser.instrument, IntegralWeakParser.topic
        )
        return info


class IntegralSPIACSParser(AlertParser):
    topic = "gcn.classic.voevent.INTEGRAL_SPIACS"
    instrument = INTEGRAL_SPI_ACS

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_integral_voevent(
            alert_msg, IntegralSPIACSParser.instrument, IntegralSPIACSParser.topic
        )
        return info
    

class IntegralWakeUpParser(AlertParser):
    topic = "gcn.classic.voevent.INTEGRAL_WAKEUP"
    instrument = INTEGRAL_MANY

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_integral_voevent(
            alert_msg, IntegralWakeUpParser.instrument, IntegralWakeUpParser.topic
        )
        return info


class IntegralRefinedParser(AlertParser):
    topic = "gcn.classic.voevent.INTEGRAL_REFINED"
    instrument = INTEGRAL_MANY

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_integral_voevent(
            alert_msg, IntegralRefinedParser.instrument, IntegralRefinedParser.topic
        )
        return info
    

class IntegralOfflineParser(AlertParser):
    topic = "gcn.classic.voevent.INTEGRAL_OFFLINE"
    instrument = INTEGRAL_MANY

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_integral_voevent(
            alert_msg, IntegralOfflineParser.instrument, IntegralOfflineParser.topic
        )
        return info
    