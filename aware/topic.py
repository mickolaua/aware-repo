__all__ = ["TOPICS", "full_topic_name_to_short", "short_topic_name_to_full"]


def full_topic_name_to_short(full_name: str) -> str:
    return full_name.lstrip("gcn.classic.voevent.")


def short_topic_name_to_full(short_name: str) -> str:
    return f"gcn.classic.voevent.{short_name}"


TOPICS = [
    "gcn.classic.voevent.FERMI_GBM_SUBTHRESH",
    "gcn.classic.voevent.SWIFT_BAT_SCALEDMAP",
    "gcn.classic.voevent.MAXI_TEST",
    "gcn.classic.voevent.MAXI_UNKNOWN",
    "gcn.classic.voevent.SWIFT_XRT_POSITION",
    "gcn.classic.voevent.SWIFT_BAT_GRB_POS_TEST",
    "gcn.classic.voevent.FERMI_GBM_FIN_INTERNAL",
    "gcn.classic.voevent.SWIFT_XRT_CENTROID",
    "gcn.classic.voevent.SWIFT_XRT_IMAGE_PROC",
    "gcn.classic.voevent.FERMI_GBM_FLT_POS",
    "gcn.classic.voevent.GRB_CNTRPART",
    "gcn.classic.voevent.ICECUBE_ASTROTRACK_GOLD",
    "gcn.classic.voevent.FERMI_LAT_MONITOR",
    "gcn.classic.voevent.FERMI_LAT_OFFLINE",
    "gcn.classic.voevent.LVC_RETRACTION",
    "gcn.classic.voevent.LVC_UPDATE",
    "gcn.classic.voevent.SWIFT_UVOT_POS_NACK",
    "gcn.classic.voevent.SNEWS",
    "gcn.classic.voevent.INTEGRAL_WEAK",
    "gcn.classic.voevent.FERMI_POINTDIR",
    "gcn.classic.voevent.SWIFT_UVOT_FCHART_PROC",
    "gcn.classic.voevent.FERMI_GBM_ALERT_INTERNAL",
    "gcn.classic.voevent.AGILE_GRB_POS_TEST",
    "gcn.classic.voevent.SWIFT_BAT_GRB_POS_ACK",
    "gcn.classic.voevent.INTEGRAL_POINTDIR",
    "gcn.classic.voevent.FERMI_GBM_POS_TEST",
    "gcn.classic.voevent.SWIFT_UVOT_DBURST_PROC",
    "gcn.classic.voevent.LVC_INITIAL",
    "gcn.classic.voevent.SWIFT_XRT_LC",
    "gcn.classic.voevent.SWIFT_XRT_IMAGE",
    "gcn.classic.voevent.SWIFT_UVOT_POS",
    "gcn.classic.voevent.SWIFT_UVOT_EMERGENCY",
    "gcn.classic.voevent.SWIFT_TOO_FOM",
    "gcn.classic.voevent.FERMI_GBM_FIN_POS",
    "gcn.classic.voevent.INTEGRAL_OFFLINE",
    "gcn.classic.voevent.SWIFT_XRT_SPECTRUM_PROC",
    "gcn.classic.voevent.FERMI_GBM_GND_INTERNAL",
    "gcn.classic.voevent.INTEGRAL_WAKEUP",
    "gcn.classic.voevent.SWIFT_XRT_THRESHPIX_PROC",
    "gcn.classic.voevent.IPN_RAW",
    "gcn.classic.voevent.MAXI_KNOWN",
    "gcn.classic.voevent.ICECUBE_ASTROTRACK_BRONZE",
    "gcn.classic.voevent.SWIFT_BAT_GRB_LC",
    "gcn.classic.voevent.AGILE_MCAL_ALERT",
    "gcn.classic.voevent.ICECUBE_CASCADE",
    "gcn.classic.voevent.SWIFT_UVOT_FCHART",
    "gcn.classic.voevent.SWIFT_XRT_SPECTRUM",
    "gcn.classic.voevent.LVC_EARLY_WARNING",
    "gcn.classic.voevent.LVC_PRELIMINARY",
    "gcn.classic.voevent.FERMI_GBM_FLT_INTERNAL",
    "gcn.classic.voevent.COINCIDENCE",
    "gcn.classic.voevent.SWIFT_BAT_QL_POS",
    "gcn.classic.voevent.SWIFT_XRT_THRESHPIX",
    "gcn.classic.voevent.SK_SN",
    "gcn.classic.voevent.FERMI_LAT_POS_TEST",
    "gcn.classic.voevent.SWIFT_POINTDIR",
    "gcn.classic.voevent.SWIFT_XRT_SPER_PROC",
    "gcn.classic.voevent.FERMI_GBM_GND_POS",
    "gcn.classic.voevent.INTEGRAL_REFINED",
    "gcn.classic.voevent.TEST_COORDS",
    "gcn.classic.voevent.SWIFT_TOO_SC_SLEW",
    "gcn.classic.voevent.KONUS_LC",
    "gcn.classic.voevent.SWIFT_ACTUAL_POINTDIR",
    "gcn.classic.voevent.HAWC_BURST_MONITOR",
    "gcn.classic.voevent.HETE_TEST",
    "gcn.classic.voevent.INTEGRAL_SPIACS",
    "gcn.classic.voevent.SWIFT_UVOT_DBURST",
    "gcn.classic.voevent.SWIFT_SC_SLEW",
    "gcn.classic.voevent.FERMI_GBM_ALERT",
    "gcn.classic.voevent.SWIFT_FOM_OBS",
    "gcn.classic.voevent.UNKNOWN",
    "gcn.classic.voevent.SWIFT_XRT_SPER",
    "gcn.classic.voevent.CALET_GBM_FLT_LC",
    "gcn.classic.voevent.SWIFT_BAT_TRANS",
    "gcn.classic.voevent.AMON_NU_EM_COINC",
    "gcn.classic.voevent.GECAM_FLT",
    "gcn.classic.voevent.GECAM_GND",
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

SUPPORTED_TOPICS = (
    LVC_TOPICS.union(INTEGRAL_TOPICS)
    .union(ICECUBE_TOPICS)
    .union(FERMI_TOPICS)
    .union(SWIFT_TOPICS)
    .union(GECAM_TOPICS)
)

TOPICS_BY_INSTRUMENT = {
    "LVC": LVC_TOPICS,
    "INTEGRAL": INTEGRAL_TOPICS,
    "ICECUBE": ICECUBE_TOPICS,
    "FERMI": FERMI_TOPICS,
    "SWIFT": SWIFT_TOPICS,
    "GECAM": GECAM_TOPICS,
}
