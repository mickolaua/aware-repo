"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
icecube.py (c) 2023
Desc: parsers of the alert messages by IceCube
Created:  2023-07-10
Modified: !date!
"""
from __future__ import annotations

from io import BytesIO, StringIO
from textwrap import dedent
from typing import Any

from astropy import units as u
from astropy.io import fits
from astropy.time import Time

from ...localization import CircularSkyMap
from ...logger import log
from ...voevent import VOEvent
from ..parser import AlertParser
from ..target_info import TargetInfo


__all__ = ["ICECUBE_ASTROTRACK_BRONZE_Parser"]


def parse_icecube_astrotrack_voevent(
    alert_msg: str | BytesIO | StringIO, instrument: str = "", topic: str = ""
) -> TargetInfo | None:
    # VOEvent document tree
    voevent = VOEvent.from_string(alert_msg)
    root = VOEvent.from_string(alert_msg).root
    log.debug(
        "Received VOEvent from %s of type %s with body \n%s",
        instrument,
        topic,
        alert_msg,
    )

    # If it GRB or not
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
    target = voevent.get_parameter_value("event_id")

    # Probability of an astrophysical origin
    P_astro = voevent.get_parameter_value("signalness", type_=float)

    # False Alarm Rate (FAR) in events per year
    FAR = voevent.get_parameter_value("FAR", type_=float)

    # Neutrino energy in TeV
    E_Nu = voevent.get_parameter_value("energy", type_=float)

    # Trigger significance
    try:
        importance = 1.0 - float(root.Why.Inference.attrib["probability"])
    except AttributeError:
        importance = 1.0

    # Collect the meta
    meta = {
        "FAR": FAR,
        "retraction": retraction,
        "signalness": P_astro,
        "event_id": target,
        "energy": E_Nu,

    }

    # Create the skymap
    if not error_radius and not ra and not dec:
        sky_map = None
    else:
        sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)

    descr = sky_map.describe() if sky_map is not None else ""
    if retraction:
        descr += "\nEVENT IS NOT A GRB"
    else:
        if FAR is not None:
            descr += f"\nFAR: {FAR:.2g} yr-1"

        if P_astro is not None:
            descr += f"\nP_Astro: {P_astro:.2g}"

        if E_Nu is not None:
            descr += f"\nE_Nu: {E_Nu:.2g} TeV"

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


class ICECUBE_ASTROTRACK_BRONZE_Parser(AlertParser):
    topic = "gcn.classic.voevent.ICECUBE_ASTROTRACK_BRONZE"
    instrument = "ICECUBE"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_icecube_astrotrack_voevent(
            alert_msg,
            ICECUBE_ASTROTRACK_BRONZE_Parser.instrument,
            ICECUBE_ASTROTRACK_BRONZE_Parser.topic,
        )
        return info


class ICECUBE_ASTROTRACK_GOLD_Parser(AlertParser):
    topic = "gcn.classic.voevent.ICECUBE_ASTROTRACK_GOLD"
    instrument = "ICECUBE"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_icecube_astrotrack_voevent(
            alert_msg,
            ICECUBE_ASTROTRACK_GOLD_Parser.instrument,
            ICECUBE_ASTROTRACK_GOLD_Parser.topic,
        )
        return info


class ICECUBE_CASCADE_Parser(AlertParser):
    topic = "gcn.classic.voevent.ICECUBE_CASCADE"
    instrument = "ICECUBE"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_icecube_astrotrack_voevent(
            alert_msg, ICECUBE_CASCADE_Parser.instrument, ICECUBE_CASCADE_Parser.topic
        )

        return info


class ICECUBE_HAWC_Parser(AlertParser):
    topic = "gcn.classic.voevent.HAWC_BURST_MONITOR"
    instrument = "ICECUBE (HAWC Burst Monitor)"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        event = VOEvent.from_string(alert_msg)
        root = VOEvent.from_string(alert_msg).root
        log.debug(
            "Received VOEvent from %s of type %s with body \n%s",
            ICECUBE_HAWC_Parser.instrument,
            ICECUBE_HAWC_Parser.topic,
            alert_msg,
        )

        info = TargetInfo(None)

        return info


class ICECUBE_NU_EM_COINC_Parser(AlertParser):
    topic = "gcn.classic.voevent.AMON_NU_EM_COINC"
    instrument = "ICECUBE"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        # VOEvent document tree
        event = VOEvent.from_string(alert_msg)
        root = VOEvent.from_string(alert_msg).root
        log.debug(
            "Received VOEvent from %s of type %s with body \n%s",
            ICECUBE_NU_EM_COINC_Parser.instrument,
            ICECUBE_NU_EM_COINC_Parser.topic,
            alert_msg,
        )

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
            isot = Time(
                astro_coords.Time.TimeInstant.ISOTime.text, format="isot"
            ).datetime
        except AttributeError:
            isot = None

        # Trigger ID
        target = event.get_parameter_value("Event_ID")

        # False Alarm Rate (FAR) in events per year
        FAR = event.get_parameter_value("FAR", type_=float)

        # Indicates that event is no longer considered astrophysical
        retraction = event.get_parameter_value("retraction", type_=bool)

        # Time search window in sec
        dt = event.get_parameter_value("deltaT", type_=float)

        # Trigger significance
        try:
            importance = 1.0 - float(root.Why.Inference.attrib["probability"])
        except AttributeError:
            importance = 1.0

        # Skymap URL
        skymap_url = event.get_parameter_value("skymap_png")

        # Collect the meta
        meta = {
            "skymap_png": skymap_url,
            "FAR": FAR,
            "retraction": retraction,
            "deltaT": dt,
            "event_id": target
        }

        # Create the skymap
        if not error_radius and not ra and not dec:
            sky_map = None
        else:
            sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)

        descr = sky_map.describe() if sky_map is not None else ""
        if retraction:
            descr += "\nEVENT IS NOT A GRB"
        else:
            if FAR is not None:
                descr += f"\nFAR: {FAR:.2g} yr-1"

            if dt is not None:
                descr += f"deltaT: {dt:.2g} sec"

            if skymap_url is not None:
                descr += f"Skymap: {skymap_url}"

        info = TargetInfo(
            sky_map,
            packet_type=ICECUBE_NU_EM_COINC_Parser.topic,
            event=target,
            origin=ICECUBE_NU_EM_COINC_Parser.instrument,
            trigger_date=isot,
            importance=importance,
            rejected=retraction,
            description=descr,
            meta=meta
        )

        return info
