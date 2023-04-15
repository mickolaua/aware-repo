from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO, StringIO
from typing import Any, Optional, Protocol

import pendulum
import requests
from astropy import units as u
from astropy.io import fits

from ..parser import AlertParser
from ..target_info import TargetInfo

from ...logger import log
from ...voevent import VOEvent
from ...localization import CircularSkyMap

__all__ = ["GBMAlertParser"]


class GBMAlertParser(AlertParser):
    topic = "gcn.classic.voevent.FERMI_GBM_FIN_POS"
    instrument = "Fermi GBM"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:

        # VOEvent document tree
        root = VOEvent.from_string(alert_msg).root
        log.debug(
            "Received VOEvent from %s of type %s with body \n%s",
            GBMAlertParser.instrument,
            GBMAlertParser.topic,
            alert_msg,
        )

        # If it GRB or not
        notgrb = root.find(".//Param[@name='Def_NOT_a_GRB']").attrib["value"]
        if notgrb == "true":
            return None

        # Astronomical coordinates and time
        obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
        astro_coords = obs_loc.AstroCoords

        # Position of the center of the localization area
        pos2d = astro_coords.Position2D
        ra = float(pos2d.Value2.C1)
        dec = float(pos2d.Value2.C2)

        # Position uncertainty (radius): stat + sys (APPROXIMATE!)
        error_radius = float(pos2d.Error2Radius) + 2.0

        # Trigger time
        isot = pendulum.parse(astro_coords.Time.TimeInstant.ISOTime.text)

        # Trigger ID
        target = root.find(".//Param[@name='TrigID']").attrib["value"]

        # Trigger significance
        importance = float(root.Why.attrib["importance"])

        # Pack target info
        sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)
        info = TargetInfo(
            sky_map,
            packet_type=GBMAlertParser.topic,
            event=target,
            origin=GBMAlertParser.instrument,
            trigger_date=isot,
            importance=importance,
        )

        return info