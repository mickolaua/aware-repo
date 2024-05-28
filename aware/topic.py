from __future__ import annotations
from typing import Any, Literal
import string

__all__ = ["TOPICS", "full_topic_name_to_short", "short_topic_name_to_full"]


VOEVENT_TOPICS = (
    "gcn.classic.voevent.AGILE_GRB_POS_TEST",
    "gcn.classic.voevent.AGILE_MCAL_ALERT",
    "gcn.classic.voevent.AMON_NU_EM_COINC",
    "gcn.classic.voevent.CALET_GBM_FLT_LC",
    "gcn.classic.voevent.COINCIDENCE",
    "gcn.classic.voevent.FERMI_GBM_ALERT_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_ALERT",
    "gcn.classic.voevent.FERMI_GBM_FIN_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_FIN_POS",
    "gcn.classic.voevent.FERMI_GBM_FLT_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_FLT_POS",
    "gcn.classic.voevent.FERMI_GBM_GND_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_GND_POS",
    "gcn.classic.voevent.FERMI_GBM_POS_TEST",
    "gcn.classic.voevent.FERMI_GBM_SUBTHRESH",
    "gcn.classic.voevent.FERMI_LAT_MONITOR",
    "gcn.classic.voevent.FERMI_LAT_OFFLINE",
    "gcn.classic.voevent.FERMI_LAT_POS_TEST",
    "gcn.classic.voevent.FERMI_POINTDIR",
    "gcn.classic.voevent.GECAM_FLT",
    "gcn.classic.voevent.GECAM_GND",
    "gcn.classic.voevent.GRB_CNTRPART",
    "gcn.classic.voevent.HAWC_BURST_MONITOR",
    "gcn.classic.voevent.HETE_TEST",
    "gcn.classic.voevent.ICECUBE_ASTROTRACK_BRONZE",
    "gcn.classic.voevent.ICECUBE_ASTROTRACK_GOLD",
    "gcn.classic.voevent.ICECUBE_CASCADE",
    "gcn.classic.voevent.INTEGRAL_OFFLINE",
    "gcn.classic.voevent.INTEGRAL_POINTDIR",
    "gcn.classic.voevent.INTEGRAL_REFINED",
    "gcn.classic.voevent.INTEGRAL_SPIACS",
    "gcn.classic.voevent.INTEGRAL_WAKEUP",
    "gcn.classic.voevent.INTEGRAL_WEAK",
    "gcn.classic.voevent.IPN_RAW",
    "gcn.classic.voevent.KONUS_LC",
    "gcn.classic.voevent.LVC_EARLY_WARNING",
    "gcn.classic.voevent.LVC_INITIAL",
    "gcn.classic.voevent.LVC_PRELIMINARY",
    "gcn.classic.voevent.LVC_RETRACTION",
    "gcn.classic.voevent.LVC_UPDATE",
    "gcn.classic.voevent.MAXI_KNOWN",
    "gcn.classic.voevent.MAXI_TEST",
    "gcn.classic.voevent.MAXI_UNKNOWN",
    "gcn.classic.voevent.SK_SN",
    "gcn.classic.voevent.SNEWS",
    "gcn.classic.voevent.SWIFT_ACTUAL_POINTDIR",
    "gcn.classic.voevent.SWIFT_BAT_GRB_LC",
    "gcn.classic.voevent.SWIFT_BAT_GRB_POS_ACK",
    "gcn.classic.voevent.SWIFT_BAT_GRB_POS_TEST",
    "gcn.classic.voevent.SWIFT_BAT_QL_POS",
    "gcn.classic.voevent.SWIFT_BAT_SCALEDMAP",
    "gcn.classic.voevent.SWIFT_BAT_TRANS",
    "gcn.classic.voevent.SWIFT_FOM_OBS",
    "gcn.classic.voevent.SWIFT_POINTDIR",
    "gcn.classic.voevent.SWIFT_SC_SLEW",
    "gcn.classic.voevent.SWIFT_TOO_FOM",
    "gcn.classic.voevent.SWIFT_TOO_SC_SLEW",
    "gcn.classic.voevent.SWIFT_UVOT_DBURST_PROC",
    "gcn.classic.voevent.SWIFT_UVOT_DBURST",
    "gcn.classic.voevent.SWIFT_UVOT_EMERGENCY",
    "gcn.classic.voevent.SWIFT_UVOT_FCHART_PROC",
    "gcn.classic.voevent.SWIFT_UVOT_FCHART",
    "gcn.classic.voevent.SWIFT_UVOT_POS_NACK",
    "gcn.classic.voevent.SWIFT_UVOT_POS",
    "gcn.classic.voevent.SWIFT_XRT_CENTROID",
    "gcn.classic.voevent.SWIFT_XRT_IMAGE_PROC",
    "gcn.classic.voevent.SWIFT_XRT_IMAGE",
    "gcn.classic.voevent.SWIFT_XRT_LC",
    "gcn.classic.voevent.SWIFT_XRT_POSITION",
    "gcn.classic.voevent.SWIFT_XRT_SPECTRUM_PROC",
    "gcn.classic.voevent.SWIFT_XRT_SPECTRUM",
    "gcn.classic.voevent.SWIFT_XRT_SPER_PROC",
    "gcn.classic.voevent.SWIFT_XRT_SPER",
    "gcn.classic.voevent.SWIFT_XRT_THRESHPIX_PROC",
    "gcn.classic.voevent.SWIFT_XRT_THRESHPIX",
    "gcn.classic.voevent.TEST_COORDS",
    "gcn.classic.voevent.UNKNOWN",
)


def alert_format_to_prefix(fmt: Literal["voevent", "binary", "json"]) -> str:
    if fmt == "voevent":
        prefix = "gcn.classic.voevent."
    elif fmt == "binary":
        prefix = "gcn.classic.binary."
    elif fmt == "json":
        prefix = "gcn.notices."
    else:
        raise ValueError(f"Unknown format `{fmt}`")

    return prefix


def full_topic_name_to_short(
    full_name: str, fmt: Literal["voevent", "binary", "json"] = "voevent"
) -> str:
    """
    Strips the prefix from the topic name.
    E.g. gcn.classic.voevent.SWIFT_XRT_SPER -> SWIFT_XRT_SPER
    """
    return full_name.lstrip(alert_format_to_prefix(fmt))


def short_topic_name_to_full(short_name: str) -> str:
    return f"gcn.classic.voevent.{short_name}"


class GCNTOpic:
    """
    A GCN alert topic

    Attributes
    ----------
    name: str
        a topic name (represents format of the message and sender)
    short_name: str
        a short name of the topic with removed prefix
        E.g. gcn.classic.voevent.SWIFT_XRT_SPER -> SWIFT_XRT_SPER
        Useful for displaying somewhere in messages.
    """

    name: str = ""
    short_name: str = ""
    __topics__: dict[str, "GCNTOpic"] = {}

    def __init__(self, **kwargs) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}(name={self.name})>"

    def __str__(self) -> str:
        return self.name

    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)
        if hasattr(cls, "name") and cls.name:
            cls.__topics__[cls.name] = cls

    @classmethod
    def get_topic(cls, name: str) -> GCNTOpic | None:
        topic = cls.__topics__.get(name, None)
        return topic
    
    @classmethod
    def list_topics(cls) -> list[GCNTOpic]:
        return [topic for (name, topic) in cls.__topics__.items() if name]

class VoeventGCNTOpic(GCNTOpic):

    @property
    def short_name(self) -> str:
        return self.name.lstrip("gcn.classic.voevent.")


class LvkJsonGCNTopic(GCNTOpic):
    """
    Common topic for all LVK alerts. Actually, duplicates several VOEvent topics.
    However, maybe more favorably than VOEvents in future.
    """

    name = "igwn.gwalert"
    short_name = "GWALERT"


class SwiftGuanoGCNTopic(GCNTOpic):
    name = "gcn.notices.swift.bat.guano"
    short_name = "SWIFT_GUANO"


class EpWXTGCNTopic(GCNTOpic):
    name = "gcn.notices.einstein_probe.wxt.alert"
    short_name = "WXT_ALERT"


class IcecubeNuemLvkSearchGCNTopic(GCNTOpic):
    name = "gcn.notices.icecube.lvk_nu_track_search"
    short_name = "LVK_NU_TRACK_SEARCH"



def safe_topic_name(topic_name) -> str:
    """
    Adjust topic name so a GCNTopic instance can be named with it. This function 
    returns the Python-safe name for an object.
    """
    # Only letters, numbers and an underscore are allowed (however, unicode characters 
    # are also valid, here we ignore them)
    safe_chars = set(string.ascii_letters + string.digits + "_")
    chars = list(topic_name)
    for i, char in enumerate(chars):
        if char not in safe_chars:
            chars[i] = "_"
    
    # Python object names can not be started with a number
    if chars[0] in set(string.digits):
        chars.insert(0, "-")

    return "".join(chars)


def create_voevent_topics() -> list[str]:
    """
    Dynamically create a list of VOEvent topics from list of strings.
    """
    voevent_topics = []
    for name in VOEVENT_TOPICS:
        tname = safe_topic_name(name.split(".")[-1])
        instance_name = tname.title() + "GCNTopic"
        topic = type(
            instance_name,
            (VoeventGCNTOpic,),
            {"name": name, "short_name": full_topic_name_to_short(name, "voevent")},
        )
        voevent_topics.append(topic())
    return voevent_topics


# Register the VOEvents
create_voevent_topics()


def get_topics() -> list[GCNTOpic]:
    """
    Get list of supported alert topics.
    """
    topics = [topic() for topic in GCNTOpic.__topics__.values()]
    return topics


TOPICS = [
    "gcn.classic.voevent.AGILE_GRB_POS_TEST",
    "gcn.classic.voevent.AGILE_MCAL_ALERT",
    "gcn.classic.voevent.AMON_NU_EM_COINC",
    "gcn.classic.voevent.CALET_GBM_FLT_LC",
    "gcn.classic.voevent.COINCIDENCE",
    "gcn.classic.voevent.FERMI_GBM_ALERT_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_ALERT",
    "gcn.classic.voevent.FERMI_GBM_FIN_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_FIN_POS",
    "gcn.classic.voevent.FERMI_GBM_FLT_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_FLT_POS",
    "gcn.classic.voevent.FERMI_GBM_GND_INTERNAL",
    "gcn.classic.voevent.FERMI_GBM_GND_POS",
    "gcn.classic.voevent.FERMI_GBM_POS_TEST",
    "gcn.classic.voevent.FERMI_GBM_SUBTHRESH",
    "gcn.classic.voevent.FERMI_LAT_MONITOR",
    "gcn.classic.voevent.FERMI_LAT_OFFLINE",
    "gcn.classic.voevent.FERMI_LAT_POS_TEST",
    "gcn.classic.voevent.FERMI_POINTDIR",
    "gcn.classic.voevent.GECAM_FLT",
    "gcn.classic.voevent.GECAM_GND",
    "gcn.classic.voevent.GRB_CNTRPART",
    "gcn.classic.voevent.HAWC_BURST_MONITOR",
    "gcn.classic.voevent.HETE_TEST",
    "gcn.classic.voevent.ICECUBE_ASTROTRACK_BRONZE",
    "gcn.classic.voevent.ICECUBE_ASTROTRACK_GOLD",
    "gcn.classic.voevent.ICECUBE_CASCADE",
    "gcn.classic.voevent.INTEGRAL_OFFLINE",
    "gcn.classic.voevent.INTEGRAL_POINTDIR",
    "gcn.classic.voevent.INTEGRAL_REFINED",
    "gcn.classic.voevent.INTEGRAL_SPIACS",
    "gcn.classic.voevent.INTEGRAL_WAKEUP",
    "gcn.classic.voevent.INTEGRAL_WEAK",
    "gcn.classic.voevent.IPN_RAW",
    "gcn.classic.voevent.KONUS_LC",
    "gcn.classic.voevent.LVC_EARLY_WARNING",
    "gcn.classic.voevent.LVC_INITIAL",
    "gcn.classic.voevent.LVC_PRELIMINARY",
    "gcn.classic.voevent.LVC_RETRACTION",
    "gcn.classic.voevent.LVC_UPDATE",
    "gcn.classic.voevent.MAXI_KNOWN",
    "gcn.classic.voevent.MAXI_TEST",
    "gcn.classic.voevent.MAXI_UNKNOWN",
    "gcn.classic.voevent.SK_SN",
    "gcn.classic.voevent.SNEWS",
    "gcn.classic.voevent.SWIFT_ACTUAL_POINTDIR",
    "gcn.classic.voevent.SWIFT_BAT_GRB_LC",
    "gcn.classic.voevent.SWIFT_BAT_GRB_POS_ACK",
    "gcn.classic.voevent.SWIFT_BAT_GRB_POS_TEST",
    "gcn.classic.voevent.SWIFT_BAT_QL_POS",
    "gcn.classic.voevent.SWIFT_BAT_SCALEDMAP",
    "gcn.classic.voevent.SWIFT_BAT_TRANS",
    "gcn.classic.voevent.SWIFT_FOM_OBS",
    "gcn.classic.voevent.SWIFT_POINTDIR",
    "gcn.classic.voevent.SWIFT_SC_SLEW",
    "gcn.classic.voevent.SWIFT_TOO_FOM",
    "gcn.classic.voevent.SWIFT_TOO_SC_SLEW",
    "gcn.classic.voevent.SWIFT_UVOT_DBURST_PROC",
    "gcn.classic.voevent.SWIFT_UVOT_DBURST",
    "gcn.classic.voevent.SWIFT_UVOT_EMERGENCY",
    "gcn.classic.voevent.SWIFT_UVOT_FCHART_PROC",
    "gcn.classic.voevent.SWIFT_UVOT_FCHART",
    "gcn.classic.voevent.SWIFT_UVOT_POS_NACK",
    "gcn.classic.voevent.SWIFT_UVOT_POS",
    "gcn.classic.voevent.SWIFT_XRT_CENTROID",
    "gcn.classic.voevent.SWIFT_XRT_IMAGE_PROC",
    "gcn.classic.voevent.SWIFT_XRT_IMAGE",
    "gcn.classic.voevent.SWIFT_XRT_LC",
    "gcn.classic.voevent.SWIFT_XRT_POSITION",
    "gcn.classic.voevent.SWIFT_XRT_SPECTRUM_PROC",
    "gcn.classic.voevent.SWIFT_XRT_SPECTRUM",
    "gcn.classic.voevent.SWIFT_XRT_SPER_PROC",
    "gcn.classic.voevent.SWIFT_XRT_SPER",
    "gcn.classic.voevent.SWIFT_XRT_THRESHPIX_PROC",
    "gcn.classic.voevent.SWIFT_XRT_THRESHPIX",
    "gcn.classic.voevent.TEST_COORDS",
    "gcn.classic.voevent.UNKNOWN",
    "gcn.notices.einstein_probe.wxt.alert",
]


LVC_TOPICS = {t for t in TOPICS if full_topic_name_to_short(t).startswith("LVC")}
INTEGRAL_TOPICS = {
    t for t in TOPICS if full_topic_name_to_short(t).startswith("INTEGRAL")
}
ICECUBE_TOPICS = {
    t for t in TOPICS if full_topic_name_to_short(t).startswith("ICECUBE")
}
FERMI_TOPICS = {t for t in TOPICS if full_topic_name_to_short(t).startswith("FERMI")}
SWIFT_TOPICS = {t for t in TOPICS if full_topic_name_to_short(t).startswith("SWIFT")}
GECAM_TOPICS = {t for t in TOPICS if full_topic_name_to_short(t).startswith("GECAM")}
EP_TOPICS = {t for t in TOPICS if full_topic_name_to_short(t).startswith("EP")}

SUPPORTED_TOPICS = (
    LVC_TOPICS.union(INTEGRAL_TOPICS)
    .union(ICECUBE_TOPICS)
    .union(FERMI_TOPICS)
    .union(SWIFT_TOPICS)
    .union(GECAM_TOPICS)
    .union(EP_TOPICS)
)

TOPICS_BY_INSTRUMENT = {
    "LVC": LVC_TOPICS,
    "INTEGRAL": INTEGRAL_TOPICS,
    "ICECUBE": ICECUBE_TOPICS,
    "FERMI": FERMI_TOPICS,
    "SWIFT": SWIFT_TOPICS,
    "GECAM": GECAM_TOPICS,
    "EP": EP_TOPICS,
}
