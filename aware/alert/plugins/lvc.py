from __future__ import annotations

from io import BytesIO, StringIO
import os
import posixpath
from textwrap import dedent
from typing import Any
import uuid

from astropy import units as u
from astropy.time import Time
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.moc import rasterize
from pathvalidate import sanitize_filename

# from astropy.utils.data import download_file

from ...localization import LVCSkyMap
from ...localization.main import lvk_uncert_level
from ...logger import log
from ...voevent import VOEvent
from ..parser import AlertParser
from ..target_info import TargetInfo
from ...config import CfgOption, dev
from ..util import is_retracted
from ...util import download_file
from ...data import cache_dir
from ...topic import full_topic_name_to_short

__all__ = ["LVC_PRELIMINARY_Parser"]


has_ns_prob_thresh = CfgOption("has_ns_prob_thresh", 0.667, float)
bbh_skymap_max_area = CfgOption("bbh_skymap_max_area", 30, float)
far_threshold_global = CfgOption("far_threshold_global", 10, float)


HZ_to_reciprocal_yr = 31557600


def parse_lvc_alert(
    alert_msg: str | BytesIO | StringIO, topic: str, create_loc: bool = True
) -> TargetInfo | None:
    rejected = True if topic == "gcn.classic.voevent.LVC_RETRACTION" else False

    root = VOEvent.from_string(alert_msg).root
    voevent = VOEvent.from_string(alert_msg)

    event_id = voevent.get_parameter_value("GraceID")
    event_page = voevent.get_parameter_value("EventPage")
    packet_type = topic
    if topic != "gcn.classic.voevent.LVC_RETRACTION":
        far = voevent.get_parameter_value("FAR", float)
        hw_inject = voevent.get_parameter_value("HardwareInj", int)
        instruments = voevent.get_parameter_value("Instruments")
        healpix_url = voevent.get_parameter_value("skymap_fits")
        nsbh_prob = voevent.get_parameter_value("NSBH", float)
        bns_prob = voevent.get_parameter_value("BNS", float)
        bbh_prob = voevent.get_parameter_value("BBH", float)
        terr_prob = voevent.get_parameter_value("Terrestrial", float)
        has_ns_prob = voevent.get_parameter_value("HasNS", float)
        has_remnant_prob = voevent.get_parameter_value("HasRemnant", float)
        pipeline = voevent.get_parameter_value("Pipeline")

    if event_id.startswith("MS"):
        if not dev.value:
            return None

    coords = voevent.root.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords
    trigger_date = Time(
        coords.Time.TimeInstant.ISOTime.text,
        format="isot",
    ).datetime

    # HEALPIX_FNAME = healpix_url.split("/")[-1].replace(",", ".")
    # resp = requests.request("GET", healpix_url)
    # healpix = BytesIO(resp.content)

    timeout = 15 * 60  # 15 min
    loc = None
    meta = {}
    try:
        if create_loc and topic != "gcn.classic.voevent.LVC_RETRACTION":
            unique = full_topic_name_to_short(topic).lower()
            local_path = posixpath.join(cache_dir.value, f"lvc_{event_id}", unique)
            os.makedirs(local_path, exist_ok=True)
            status = download_file(healpix_url, path=local_path, timeout=timeout)
            filename = posixpath.split(healpix_url)[1]
            sky_map = read_sky_map(
                posixpath.join(local_path, filename), nest=True, moc=True
            )
            meta = sky_map.meta
            data = sky_map
            uniq = data["UNIQ"]
            probdensity = data["PROBDENSITY"]
            distmu = meta["distmean"] * u.Mpc
            distsigma = meta["diststd"] * u.Mpc

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
                event_page=event_page,
            )
    except Exception as e:
        log.error("Skymap is not available for LVK %s", event_id, exc_info=e)

    if topic != "gcn.classic.voevent.LVC_RETRACTION":
        # Do not send spurious detections
        if far * HZ_to_reciprocal_yr > far_threshold_global.value:
            if not dev.value:
                return None
        
        # We do not want to observe BBH event, which do not generate electromagnetic
        # counterparts. But let's check a small BBH localization for a possible signal.
        create_loc = create_loc and terr_prob < 0.333 and (
            (has_ns_prob > has_ns_prob_thresh.value and bbh_prob < 0.333)
            or 
            (loc.area() < bbh_skymap_max_area.value * u.deg * u.deg)
        )

    # Retraction alerts has only a few parameters compared with other ones
    if topic != "gcn.classic.voevent.LVC_RETRACTION":
        importance = has_ns_prob
        description = dedent(
            f"""
            FAR: {far:.3g} Hz (1 per {1/far/HZ_to_reciprocal_yr:.3g} yr^-1)
            P_BNS: {bns_prob:.1g}
            P_NSBH: {nsbh_prob:.3g}
            P_BBH: {bbh_prob:.3g}
            P_Terr: {terr_prob:.3g}
            P_hasNS: {has_ns_prob:.3g}
            P_hasRemnant: {has_remnant_prob:.3g}
            H/W injection: {hw_inject:d}
            GraceDB URL: {event_page:s}
            Skymap URL: {healpix_url:s}
            Instruments: {instruments:s}
            Algorithm: {pipeline:s}
            """
        )

        if loc:
            description += f"\n{lvk_uncert_level.value*100:.3g}% area: {loc.area():.3g}"
            if meta:
                description += (
                    f"\nDistance: {distmu:.3g} +/- {2*distsigma:.3g} (2sigma)"
                )
        else:
            description += "\nSkymap is not available at the moment."
    else:
        description = dedent(
            f"""
            GraceDB URL: {event_page:s}
            """
        )
        importance = 0.0

    info = TargetInfo(
        loc,
        packet_type=packet_type,
        event=event_id,
        origin="LVC",
        meta=meta,
        importance=importance,
        rejected=rejected,
        trigger_date=trigger_date,
        description=description,
    )

    return info


class LVC_EARLY_WARNING_Parser(AlertParser):
    """Parser of early warning LVK alert messages"""

    topic: str = "gcn.classic.voevent.LVC_EARLYWARNING"
    instrument: str = "LVC"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_lvc_alert(alert_msg, LVC_EARLY_WARNING_Parser.topic)


class LVC_PRELIMINARY_Parser(AlertParser):
    """Parser of preliminary LVK alert messages"""

    topic: str = "gcn.classic.voevent.LVC_PRELIMINARY"
    instrument: str = "LVC"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_lvc_alert(alert_msg, LVC_PRELIMINARY_Parser.topic)


class LVC_INITIAL_Parser(LVC_PRELIMINARY_Parser):
    """Parser of initial LVK alert messages"""

    topic: str = "gcn.classic.voevent.LVC_INITIAL"
    instrument: str = "LVC"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_lvc_alert(alert_msg, LVC_INITIAL_Parser.topic)


class LVC_RETRACTION_Parser(LVC_PRELIMINARY_Parser):
    """Parser of retracked LVK alert messages"""

    topic: str = "gcn.classic.voevent.LVC_RETRACTION"
    instrument: str = "LVC"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_lvc_alert(alert_msg, LVC_RETRACTION_Parser.topic, create_loc=False)


class LVC_UPDATE_Parser(LVC_PRELIMINARY_Parser):
    """Parser of update LVK alert messages"""

    topic: str = "gcn.classic.voevent.LVC_UPDATE"
    instrument: str = "LVC"

    @staticmethod
    def parse_alert(
        alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
    ) -> TargetInfo | None:
        return parse_lvc_alert(alert_msg, LVC_UPDATE_Parser.topic)
