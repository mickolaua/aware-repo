from __future__ import annotations
from functools import partial

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
from aware.util import render_number_as_rich_utf8


__all__ = ["MAXIKnownAlertParser"]


MAXI = "MAXI"


def parse_maxi_alert(
    alert_msg: str | BytesIO | StringIO, instr: str, topic: str
) -> TargetInfo:
    # VOEvent document tree
    voevent = VOEvent.from_string(alert_msg)
    root = voevent.root
    log.debug(
        "Received VOEvent from %s of type %s with body \n%s",
        instr,
        topic,
        alert_msg,
    )

    # If it GRB or not
    notgrb = voevent.get_parameter_value("Def_NOT_a_GRB", type_=bool)
    if notgrb:
        return None

    # Astronomical coordinates and time
    obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
    astro_coords = obs_loc.AstroCoords

    # Position of the center of the localization area
    pos2d = astro_coords.Position2D
    ra = float(pos2d.Value2.C1)
    dec = float(pos2d.Value2.C2)

    # Position uncertainty (radius): stat only
    error_radius = float(pos2d.Error2Radius)

    # Trigger time
    isot = Time(astro_coords.Time.TimeInstant.ISOTime.text, format="isot").datetime

    # Trigger ID
    target = root.find(".//Param[@name='TrigID']").attrib["value"]

    # Trigger significance
    importance = float(root.Why.attrib["importance"])

    # Pack target info
    sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)

    # Meta data
    meta = {
        "Source_Name": voevent.get_parameter_value("Source_Name"),
        "Source_Flux": voevent.get_parameter_value("Source_Flux", type_=float) # mCrab
    }

    # Description
    descr = ""
    if name := meta["Source_Name"]:
        descr += f"\nName: {name}\n"

    if flux := meta["Source_Flux"]:
        flux_erg = flux * 2.4e-11
        flux_erg_str = render_number_as_rich_utf8(flux_erg, precision=3)
        descr += f"Flux: {flux:.3f} mCrab ({flux_erg_str} erg cm\u207B\u00B2 s\u207B\u00B9)"

    info = TargetInfo(
        sky_map,
        packet_type=topic,
        event=target,
        origin=instr,
        trigger_date=isot,
        importance=importance,
        description=descr,
        meta=meta
    )

    return info


# class GBMAlertAlertParser(AlertParser):
#     topic = "gcn.classic.voevent.FERMI_GBM_ALERT"
#     instrument = FERMI_GBM

#     @staticmethod
#     def parse_alert(
#         alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
#     ) -> TargetInfo | None:
#         # VOEvent document tree
#         voevent = VOEvent.from_string(alert_msg)
#         root = voevent.root
#         log.debug(
#             "Received VOEvent from %s of type %s with body \n%s",
#             GBMAlertAlertParser.instrument,
#             GBMAlertAlertParser.topic,
#             alert_msg,
#         )
#         # Trigger time
#         isot = Time(
#             (
#                 root
#                 .WhereWhen
#                 .ObsDataLocation
#                 .ObservationLocation
#                 .AstroCoords
#                 .Time
#                 .TimeInstant
#                 .ISOTime
#                 .text
#             ),
#             format="isot",
#         ).datetime
#         importance = voevent.get_parameter_value("Importance", type_=float)
#         event = voevent.get_parameter_value("TrigID")
#         lc_url = voevent.get_parameter_value("LightCurve_URL")

#         meta = {
#             "Trig_Signif": voevent.get_parameter_value("Trig_Signif", type_=float),
#             "Trig_Dur": voevent.get_parameter_value("Trig_Dur", type_=float),

#         }

#         descr = f"""
# {GBMAlertAlertParser.instrument} has registered candidate transient event.
# Wait for upcoming messages, which could reveal the nature of the event.
# """
#         if lc_url:
#             descr += f"Light Curve URL: {lc_url}\n"

#         if T90 := meta["Trig_Dur"]:
#             descr += f"Trigger Duration: {T90:.3g} s\n"

#         if sigma := meta["Trig_Signif"]:
#             descr += f"Trigger Significance: {sigma:.3g}\n"

#         info = TargetInfo(
#             None,
#             packet_type=GBMAlertAlertParser.topic,
#             origin=GBMAlertAlertParser.instrument,
#             trigger_date=isot,
#             description=descr,
#             importance=importance,
#             event=event,
#             meta=meta
#         )

#         return info


class MAXIUnknownAlertParser(AlertParser):
    topic = "gcn.classic.voevent.MAXI_UNKNOWN"
    instrument = MAXI

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_maxi_alert(
            alert_msg,
            MAXIUnknownAlertParser.instrument,
            MAXIUnknownAlertParser.topic,
        )
        return info
    

class MAXIKnownAlertParser(AlertParser):
    topic = "gcn.classic.voevent.MAXI_KNOWN"
    instrument = MAXI

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        info = parse_maxi_alert(
            alert_msg,
            MAXIKnownAlertParser.instrument,
            MAXIKnownAlertParser.topic,
        )
        return info
