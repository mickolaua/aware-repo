from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO, StringIO
from typing import Any, Optional, Protocol

import pendulum
import requests
from astropy import units as u
from astropy.io import fits
from ligo.skymap.moc import rasterize
from ligo.skymap.io import read_sky_map

from ..parser import AlertParser
from ..target_info import TargetInfo
from ...logger import log
from ...voevent import VOEvent
from ...localization import LVCSkyMap


__all__ = ["LVC_PRELIMINARY_Parser"]


class LVC_PRELIMINARY_Parser(AlertParser):
    """Parser of preliminary LVK alert messages"""

    topic: str = "gcn.classic.voevent.LVC_PRELIMINARY"
    instrument: str = "LVC"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        root = VOEvent.from_string(alert_msg).root

        event_id = "Unknown"
        packet_type = 150
        bns_prob = 0
        nsbh_prob = 0
        bbh_prob = 0
        terr_prob = 0
        has_ns_prob = 0
        has_remnant_prob = 0
        has_remnant_prob = 0
        hw_inject = 1
        healpix_url = ""
        event_page = ""

        for p in root.findall(".//What/Param"):
            if p.attrib["name"] == "GraceID":
                event_id = p.attrib["value"]

            if p.attrib["name"] == "HardwareInj":
                hw_inject = int(p.attrib["value"])

            if p.attrib["name"] == "EventPage":
                event_page = p.attrib["value"]

        for p in root.findall(".//What/Group/Param"):
            if p.attrib["name"] == "Packet_Type":
                packet_type = int(p.attrib["value"])

            if p.attrib["name"] == "BNS":
                bns_prob = float(p.attrib["value"])

            if p.attrib["name"] == "NSBH":
                nsbh_prob = float(p.attrib["value"])

            if p.attrib["name"] == "BBH":
                bbh_prob = float(p.attrib["value"])

            if p.attrib["name"] == "Terrestrial":
                terr_prob = float(p.attrib["value"])

            if p.attrib["name"] == "HasNS":
                has_ns_prob = float(p.attrib["value"])

            if p.attrib["name"] == "HasRemnant":
                has_remnant_prob = float(p.attrib["value"])

            if p.attrib["name"] == "HardwareInj":
                hw_inject = int(p.attrib["value"])

            if p.attrib["name"] == "skymap_fits":
                healpix_url = p.attrib["value"]


        trigger_date = pendulum.parse(root.find(".//ISOTime").text)

        # HEALPIX_FNAME = healpix_url.split("/")[-1].replace(",", ".")
        # resp = requests.request("GET", healpix_url)
        # healpix = BytesIO(resp.content)

        sky_map= read_sky_map(healpix_url, moc=True)
        hdr = sky_map.meta
        data = sky_map
        # with fits.open(healpix, mode="readonly") as hdul:
        #     hdr = hdul[1].header
        #     data = hdul[1].data

        uniq = data["UNIQ"]
        probdensity = data["PROBDENSITY"]
        distmu = hdr["distmean"] * u.Mpc
        distsigma = hdr["diststd"] * u.Mpc

        loc = LVCSkyMap(
            bns_prob=bns_prob,
            nsbh_prob=nsbh_prob,
            bbh_prob=bbh_prob,
            terr_prob=terr_prob,
            has_ns_prob=has_ns_prob,
            has_remnant_prob=has_remnant_prob,
            hw_inject=hw_inject,
            uniq=uniq, 
            data=data,
            probdensity=probdensity,
            distmu=distmu,
            distsigma=distsigma,
            packet_type=packet_type,
            event_page=event_page
        )

        # It is considered only BNS or NSBH can produce electromagnetic
        # transient in which we interested
        importance = has_ns_prob
        info = TargetInfo(
            loc,
            event=event_id,
            origin="LVC",
            meta={},
            importance=importance,
            trigger_date=trigger_date,
        )

        return info