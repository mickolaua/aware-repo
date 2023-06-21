from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from io import BytesIO, StringIO
from typing import Any, Optional, Protocol

import pendulum
import requests
from astropy import units as u
from astropy.io import fits

from ..logger import log
from ..voevent import VOEvent
from ..localization import CircularSkyMap, Localization, LVCSkyMap


__all__ = ["alert_parsers", "TargetInfo"]




# class GECAM_FLT_Parser(AlertParser):
#     topic: str = "gcn.classic.voevent.GECAM_FLT"
#     instrument: str = "GECAM"

#     @staticmethod
#     def parse_alert(
#         alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
#     ) -> NewTargetInfo | None:
#         info = parse_voevent(alert_msg, "GECAM FLIGHT", name_field="Trigger_Number")
#         info.origin = GECAM_FLT_Parser.instrument
#         return info


# class GECAM_GND_Parser(AlertParser):
#     topic: str = "gcn.classic.voevent.GECAM_GND"
#     instrument: str = "GECAM"

#     @staticmethod
#     def parse_alert(
#         alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
#     ) -> NewTargetInfo | None:
#         info = parse_voevent(alert_msg, "GECAM GROUND", name_field="Trigger_Number")
#         info.origin = GECAM_GND_Parser.instrument
#         return info


# class SWIFTBATAlertParser(AlertParser):
#     topic = "gcn.classic.voevent.SWIFT_BAT_GRB_POS_ACK"
#     instrument = "Swift BAT"

#     @staticmethod
#     def parse_alert(
#         alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
#     ) -> NewTargetInfo | None:

#         # VOEvent document tree
#         root = VOEvent.from_string(alert_msg).root

#         # We are only interested in GRBs
#         is_grb = root.find(".//Param[@name='GRB_Identified']").attrib["value"]
#         if is_grb == "false":
#             return None

#         # Astronomical coordinates and time
#         obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
#         astro_coords = obs_loc.AstroCoords

#         # Position of the center of the localization area
#         pos2d = astro_coords.Position2D
#         ra = float(pos2d.Value2.C1)
#         dec = float(pos2d.Value2.C2)

#         # Position uncertainty (radius)
#         error_radius = float(pos2d.Error2Radius)

#         # Trigger time
#         isot = pendulum.parse(astro_coords.Time.TimeInstant.ISOTime.text)

#         # Trigger ID
#         # target = root.find(".//Param[@name='TrigID']").attrib["value"]
#         target = root.find(".//Param[@name='TrigID']").attrib["value"]

#         # Trigger significance
#         importance = float(root.Why.attrib["importance"])

#         # Pack target info
#         sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)
#         info = NewTargetInfo(
#             sky_map,
#             packet_type=SWIFTBATAlertParser.topic,
#             event=target,
#             origin=SWIFTBATAlertParser.instrument,
#             trigger_date=isot,
#             importance=importance,
#         )

#         return info


# class SWIFTXRTAlertParser(AlertParser):
#     topic = "gcn.classic.voevent.SWIFT_XRT_POSITION"
#     instrument = "Swift XRT"

#     @staticmethod
#     def parse_alert(
#         alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
#     ) -> NewTargetInfo | None:

#         # VOEvent document tree
#         root = VOEvent.from_string(alert_msg).root

#         # Astronomical coordinates and time
#         obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
#         astro_coords = obs_loc.AstroCoords

#         # Position of the center of the localization area
#         pos2d = astro_coords.Position2D
#         ra = float(pos2d.Value2.C1)
#         dec = float(pos2d.Value2.C2)

#         # Position uncertainty (radius)
#         error_radius = float(pos2d.Error2Radius)

#         # Trigger time
#         isot = pendulum.parse(astro_coords.Time.TimeInstant.ISOTime.text)

#         # Trigger ID
#         # target = root.find(".//Param[@name='TrigID']").attrib["value"]
#         target = root.find(".//Param[@name='TrigID']").attrib["value"]

#         # Trigger significance
#         importance = float(root.Why.attrib["importance"])

#         # Pack target info
#         sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)
#         info = NewTargetInfo(
#             sky_map,
#             packet_type=SWIFTXRTAlertParser.topic,
#             event=target,
#             origin=SWIFTXRTAlertParser.instrument,
#             trigger_date=isot,
#             importance=importance,
#         )

#         return info


# class SWIFTUVOTAlertParser(AlertParser):
#     topic = "gcn.classic.voevent.SWIFT_UVOT_POS"
#     instrument = "Swift UVOT"

#     @staticmethod
#     def parse_alert(
#         alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
#     ) -> NewTargetInfo | None:

#         # VOEvent document tree
#         root = VOEvent.from_string(alert_msg).root

#         # Astronomical coordinates and time
#         obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
#         astro_coords = obs_loc.AstroCoords

#         # Position of the center of the localization area
#         pos2d = astro_coords.Position2D
#         ra = float(pos2d.Value2.C1)
#         dec = float(pos2d.Value2.C2)

#         # Position uncertainty (radius)
#         error_radius = float(pos2d.Error2Radius)

#         # Trigger time
#         isot = pendulum.parse(astro_coords.Time.TimeInstant.ISOTime.text)

#         # Trigger ID
#         # target = root.find(".//Param[@name='TrigID']").attrib["value"]
#         target = root.find(".//Param[@name='TrigID']").attrib["value"]

#         # Trigger significance
#         importance = float(root.Why.attrib["importance"])

#         # Pack target info
#         sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)
#         info = NewTargetInfo(
#             sky_map,
#             packet_type=SWIFTUVOTAlertParser.topic,
#             event=target,
#             origin=SWIFTUVOTAlertParser.instrument,
#             trigger_date=isot,
#             importance=importance,
#         )

#         return info


# class SWIFT_ACTUAL_POINTDIR_AlertParser(AlertParser):
#     topic = "gcn.classic.voevent.SWIFT_ACTUAL_POINTDIR"
#     instrument = "Swift BAT"

#     @staticmethod
#     def parse_alert(
#         alert_msg: str | BytesIO | StringIO, *args: Any, **kwargs: Any
#     ) -> NewTargetInfo | None:

#         # VOEvent document tree
#         root = VOEvent.from_string(alert_msg).root

#         # Astronomical coordinates and time
#         obs_loc = root.WhereWhen.ObsDataLocation.ObservationLocation
#         astro_coords = obs_loc.AstroCoords

#         # Position of the center of the localization area
#         pos2d = astro_coords.Position2D
#         ra = float(pos2d.Value2.C1)
#         dec = float(pos2d.Value2.C2)

#         # Position uncertainty (radius)
#         error_radius = float(pos2d.Error2Radius)

#         # Trigger time
#         isot = pendulum.parse(astro_coords.Time.TimeInstant.ISOTime.text)

#         # Trigger ID
#         # target = root.find(".//Param[@name='TrigID']").attrib["value"]
#         target = root.find(".//Param[@name='TrigID']").attrib["value"]

#         # Pack target info
#         sky_map = CircularSkyMap(ra_center=ra, dec_center=dec, radius=error_radius)
#         info = NewTargetInfo(
#             sky_map,
#             packet_type=SWIFT_ACTUAL_POINTDIR_AlertParser.topic,
#             event=target,
#             origin=SWIFT_ACTUAL_POINTDIR_AlertParser.instrument,
#             trigger_date=isot,
#             importance=0.5,
#         )

#         return info


# alert_parsers: dict[str, AlertParser] = {
#     GBMAlertParser.topic: GBMAlertParser,
#     # GECAM_FLT_Parser.topic: GECAM_FLT_Parser,
#     # GECAM_GND_Parser.topic: GECAM_GND_Parser,
#     # SWIFT_BAT_GRB_POS_TEST_Parser.topic: SWIFT_BAT_GRB_POS_TEST_Parser,
#     # SWIFT_BAT_GRB_POS_ACK_Parser.topic: SWIFT_BAT_GRB_POS_ACK_Parser,
#     # SWIFT_XRT_POSITION_Parser.topic: SWIFT_XRT_POSITION_Parser,
#     LVC_PRELIMINARY_Parser.topic: LVC_PRELIMINARY_Parser,
#     SWIFTBATAlertParser.topic: SWIFTBATAlertParser,
#     SWIFTXRTAlertParser.topic: SWIFTXRTAlertParser,
#     SWIFTUVOTAlertParser.topic: SWIFTUVOTAlertParser,
#     SWIFT_ACTUAL_POINTDIR_AlertParser.topic: SWIFT_ACTUAL_POINTDIR_AlertParser,
# }
