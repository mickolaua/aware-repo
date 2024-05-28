from __future__ import annotations

import os
import posixpath
from io import BytesIO, StringIO
from textwrap import dedent
from typing import Any

import numpy as np
from astropy import units as u
from astropy.time import Time
from ligo.skymap.io.fits import read_sky_map

import aware.data

from ...config import CfgOption, dev
from ...localization import LVCSkyMap
from ...localization.main import lvk_uncert_level
from ...logger import log
from ...topic import full_topic_name_to_short
from ...util import download_file, render_number_as_rich_utf8
from ...voevent import VOEvent
from ..parser import AlertParser
from ..target_info import TargetInfo

__all__ = [
    "LVC_PRELIMINARY_Parser",
    "LVC_EARLY_WARNING_Parser",
    "LVC_INITIAL_Parser",
    "LVC_RETRACTION_Parser",
    "LVC_UPDATE_Parser",
]


has_ns_prob_thresh = CfgOption("has_ns_prob_thresh", 0.667, float)
bbh_skymap_max_area = CfgOption("bbh_skymap_max_area", 30, float)
far_threshold_global = CfgOption("far_threshold_global", 10, float)


HZ_to_reciprocal_yr = 31557600


def is_observable_bbh(loc: LVCSkyMap) -> bool:
    """Is the localization of the BBH event should be observed? Here, we attempt to
    find a likely possible optical transients occured due to the capture of the
    surrounding material to the merger of two black holes.

    Parameters
    ----------
    loc : LVCSkyMap
        a localization skymap

    Returns
    -------
    bool
        True if observable else False
    """
    return loc.area() < bbh_skymap_max_area.value * u.deg * u.deg


def has_ns(has_ns_prob: float, bbh_prob: float) -> bool:
    """Estimate the probability of having the neutron star in the merger.

    Parameters
    ----------
    has_ns_prob : float
        a probabilty to have a neutron star, given by LVK
    bbh_prob : float
        a probabilty to classify the merger as a BBH, given by LVK

    Returns
    -------
    bool
        A probability of having the neutron star in the merger in range 0..1
    """
    return has_ns_prob > has_ns_prob_thresh.value and bbh_prob < 0.333


def is_astro(terr_prob: float) -> bool:
    """Estimate the probability to classify the merger as an astrophysical event.

    Parameters
    ----------
    terr_prob : float
        a probabilty to classify the merger as a glitch

    Returns
    -------
    bool
        A probability to classify the merger as an astrophysical event in range 0..1
    """
    return terr_prob < 0.333


def calc_importance(x: float) -> float:
    """Calculates importance of the given LVK event

    Parameters
    ----------
    x : float
        FAR in units of events per year

    Returns
    -------
    float
        importance in range from 0 to 1
    """
    # norm_decay chosen in such way that if x = 1.0 event per year, importance
    # will be about 0.99
    norm_decay = 100.0
    1.0 if x < 0.0 else 1.0 if x > 1e3 else np.exp(-x / norm_decay)


def parse_lvc_alert(
    alert_msg: str | BytesIO | StringIO, topic: str, create_loc: bool = True
) -> TargetInfo | None:
    """Parse LVC alerts.

    Parameters
    ----------
    alert_msg : str | BytesIO | StringIO
        a raw LVC alert message in the VOEvent format
    topic : str
        an alert type (topic in terms of Kafka)
    create_loc : bool, optional
        create localization skymap, by default True

    Returns
    -------
    TargetInfo | None
        The target info or None in case of a mock data event
    """

    rejected = True if topic == "gcn.classic.voevent.LVC_RETRACTION" else False

    root = VOEvent.from_string(alert_msg).root
    voevent = VOEvent.from_string(alert_msg)

    event_id = voevent.get_parameter_value("GraceID")
    event_page = voevent.get_parameter_value("EventPage")
    packet_type = topic

    far = voevent.get_parameter_value("FAR", float, default=0.0)
    hw_inject = voevent.get_parameter_value("HardwareInj", int, default=0)
    instruments = voevent.get_parameter_value("Instruments", default="")
    healpix_url = voevent.get_parameter_value("skymap_fits", default="")
    nsbh_prob = voevent.get_parameter_value("NSBH", float, default=0.0)
    bns_prob = voevent.get_parameter_value("BNS", float, default=0.0)
    bbh_prob = voevent.get_parameter_value("BBH", float, default=0.0)
    terr_prob = voevent.get_parameter_value("Terrestrial", float, default=0.0)
    has_ns_prob = voevent.get_parameter_value("HasNS", float, default=0.0)
    has_remnant_prob = voevent.get_parameter_value("HasRemnant", float, default=0.0)
    pipeline = voevent.get_parameter_value("Pipeline", default="")

    if event_id.startswith("MS"):
        if not dev.value:
            return None

    coords = voevent.root.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords
    trigger_date = Time(
        coords.Time.TimeInstant.ISOTime.text,
        format="isot",
    ).datetime

    timeout = 15 * 60  # 15 min
    loc = None
    meta = {}
    tries = 5
    for i in range(1, tries + 1):
        try:
            if create_loc and topic != "gcn.classic.voevent.LVC_RETRACTION":
                log.debug(
                    "Downloading HEALPix for LVK %s; try %d/%d", event_id, i, tries
                )
                unique = full_topic_name_to_short(topic).lower()
                local_path = os.path.join(
                    aware.data.cache_dir.value, f"lvc_{event_id}", unique
                )
                os.makedirs(local_path, exist_ok=True)
                status = download_file(healpix_url, path=local_path, timeout=timeout)
                filename = posixpath.split(healpix_url)[1]
                sky_map = read_sky_map(
                    os.path.join(local_path, filename), nest=True, moc=True
                )
                meta = sky_map.meta
                data = sky_map
                uniq = data["UNIQ"]
                probdensity = data["PROBDENSITY"]
                distmu = meta.get("distmean", 0.0) * u.Mpc
                distsigma = meta.get("diststd", 0.0) * u.Mpc

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
                if status:
                    break
        except Exception as e:
            if i == tries:
                log.error("Failed to retrieve skymap for LVK %s", event_id, exc_info=e)

    # Do not send spurious detections
    if far * HZ_to_reciprocal_yr > far_threshold_global.value:
        if not dev.value:
            return None

    # We do not want to observe BBH event, which do not generate electromagnetic
    # counterparts. But let's check a small BBH localization for a possible signal.
    create_loc = (
        create_loc
        and is_astro(terr_prob)
        and (has_ns(has_ns_prob, bbh_prob) or is_observable_bbh(loc))
    )

    # Retraction alerts has only a few parameters compared with other ones
    x = far * HZ_to_reciprocal_yr

    if topic == "gcn.classic.voevent.LVC_RETRACTION":
        importance = 0.0
    else:
        importance = calc_importance(x)

    if far:
        far_str = render_number_as_rich_utf8(far, precision=3)
        repr_far_str = render_number_as_rich_utf8(
            1 / far / HZ_to_reciprocal_yr, precision=3
        )
    else:
        far_str = "0.0"
        repr_far_str = "\u221E"
    description = dedent(
        f"""
        FAR: {far_str} Hz (1 per {repr_far_str} yr)
        P_BNS: {bns_prob:.2f}
        P_NSBH: {nsbh_prob:.2f}
        P_BBH: {bbh_prob:.2f}
        P_Terr: {terr_prob:.2f}
        P_hasNS: {has_ns_prob:.2f}
        P_hasRemnant: {has_remnant_prob:.2f}
        H/W injection: {hw_inject:d}
        GraceDB URL: {event_page:s}
        Skymap URL: {healpix_url:s}
        Instruments: {instruments:s}
        Algorithm: {pipeline:s}
        """
    )

    if meta:
        # distmu_str = render_number_as_rich_utf8(distmu.value, 3)
        # distsigma_str = render_number_as_rich_utf8(2 * distsigma.value, 3)
        description += f"Distance: {distmu.value:.1f} \u00B1 {distsigma.value:.1f} Mpc (2\u03C3)"
    else:
        description += "Skymap is not available at the moment."

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
