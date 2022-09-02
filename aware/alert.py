from __future__ import annotations
from abc import ABCMeta
import atexit
from dataclasses import dataclass
from datetime import datetime
import importlib
from inspect import isabstract
import sys
from typing import Any, Optional, Protocol
from confluent_kafka import Message
from lxml import etree
from io import BytesIO, StringIO
from astropy.time import Time

from .plugin import hookimpl, hookspec, plugin_manager
import os


__all__ = ["alert_parsers", "TargetInfo"]


# Registry for alert parser by topic
class AlertParsers:
    def __init__(self, **kwargs: dict[str, AlertParser]) -> None:
        self._parsers = {}
        self._parsers.update(kwargs)
    
    def _is_parser_exists(self, parser_id: str) -> bool:
        return parser_id in self._parsers

    def add_parser(self, parser_id: str, parser: AlertParser) -> None:
        if self._is_parser_exists(parser_id):
            raise ValueError(
                f"Alert parser with id: '{parser_id}' already exists.")
        else:
            self._parsers[parser_id] = parser

    def __setitem__(self, parser_id: str, parser: AlertParser) -> None:
        self.add_parser(parser_id, parser)

    def __getitem__(self, parser_id: str) -> None:
        self.get_parser(parser_id)

    def get_parser(self, parser_id: str) -> None:
        if self._is_parser_exists(parser_id):
            return self._parsers.get(parser_id)
        else:
            raise ValueError(
                f"Alert parser with id: '{parser_id}' is not exist.")

    def remove_parser(self, parser_id: str) -> None:
        if self._is_parser_exists(parser_id):
            self._parsers.pop(parser_id)
        else:
            raise ValueError(
                f"Alert parser with id: '{parser_id}' is not exist.")

    def __delitem__(self, parser_id: str) -> None:
        self.remove_parser(parser_id)

    @property
    def parsers(self) -> list[str]:
        return list(self._parsers.keys())

    def __str__(self):
        return str(self._parsers)

    def __repr__(self):
        return repr(self._parsers)
            

alert_parsers = AlertParsers()
    

def _get_xml_root(msg: str | BytesIO | StringIO):
    return etree.fromstring(msg)


@dataclass
class TargetInfo:
    ra_center: Optional[float] = None
    dec_center: Optional[float] = None
    error_radius1: Optional[float] = None
    error_radius2: Optional[float] = None
    event: Optional[str] = None
    origin: Optional[str] = None
    trigger_date: Optional[datetime] = None
    importance: Optional[float] = None


# @hookspec
# def parse_alert(alert_msg: str | BytesIO | StringIO, 
#                 *args, **kwargs) -> TargetInfo | None:
#     """
#     Hook for parsing alert message and retrieving the target info 
#     (i.e. coordinates, name, UTC, etc.). This is hook specification, concrete 
#     hook implementation must be implemented as a class with parse_alert method 
#     or just regular function. Both must be wrapped with hookimpl wrapper. 

#     Parameters
#     ----------
#     alert_msg : str | BytesIO | StringIO
#         an alert message

#     Returns
#     -------
#     TargetInfo
#         info on alert target (object of interest -- the cause of the alert)

#     Raises
#     ------
#     This hook **specification** raises NotImplementedError
#     """
#     raise NotImplementedError


class AlertParser(ABCMeta):
    topic: str = ""

    @staticmethod
    def parse_alert(alert_msg: str | BytesIO | StringIO, 
                    *args: Any, **kwargs: Any) -> TargetInfo | None:
        raise NotImplementedError


class FERMI_GBM_FIN_POS_AlertParser(AlertParser):
    topic: str = "gcn.classic.voevent.FERMI_GBM_FIN_POS"

    @staticmethod
    def parse_alert(alert_msg: str | BytesIO | StringIO, 
                    *args: Any, **kwargs: Any) -> TargetInfo | None:
        # VOEvent document tree
        root = _get_xml_root(alert_msg)
        
        # Check that alert were sent by Fermi GBM
        descr = root.find("./How/Description").text
        if descr != "Fermi Satellite, GBM Instrument":
            raise ValueError("This is not correct Fermi GBM VOEvent")

        # Position of the center of the localization area
        pos2d = root.find(".//{*}Position2D")
        ra = float(pos2d.find("./Value2/C1").text) 
        dec = float(pos2d.find("./Value2/C2").text)

        # Position uncertainty (radius)
        error_radius = float(pos2d.find("./Error2Radius").text)

        # GBM trigger time
        isot = Time(root.find(".//{*}ISOTime").text, format="isot").datetime

        # GBM triger ID
        target = ""
        for param in root.findall('./What/Param'):
            if param.attrib["name"] == "TrigID":
                target = param.attrib["value"]

        # Trigger significance
        importance = float(root.find("./Why").attrib["importance"])
        
        # Pack target info
        info =  TargetInfo(ra_center=ra, dec_center=dec, 
                           error_radius1=error_radius, 
                           error_radius2=error_radius,
                           event=target, importance=importance, 
                           trigger_date=isot, origin=descr)

        return info


# def registrate_all_parsers():
#     for obj in dir(__dict__):
#         try:
#             if issubclass(obj, AlertParser) and obj is not AlertParser:
#                 alert_parsers.add_parser(obj.topic, obj)
#         except TypeError:
#             pass


# registrate_all_parsers()

# atexit.register(registrate_all_parsers)
# alert_parsers.add_parser(FERMI_GBM_FIN_POS_AlertParser.topic,
#                          FERMI_GBM_FIN_POS_AlertParser)

