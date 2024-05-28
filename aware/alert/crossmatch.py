from __future__ import annotations
from typing import Sequence

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time, TimeDelta
from mocpy import MOC
from sqlalchemy import Row

from aware.localization.main import Localization

from .target_info import TargetInfo
from ..sql.models import Alert, create_session, dbconnect
from ..config import CfgOption
from ..logger import log
from ..angle import coord2str
import pickle
import aiomisc


# Maximal date difference to match a pair of alerts
max_date_diff_min = CfgOption("max_date_diff_min", 5.0, float)


def crossmatch_alerts(
    info: TargetInfo, max_date_diff: u.Unit = max_date_diff_min.get_value() * u.min
):
    log.debug("Cross-matching alert with a database ...")
    event = info.event
    date = info.trigger_date
    origin = info.origin
    moc = info.localization.moc()
    center = moc.barycenter()
    log.debug(
        "\nEvent: %s\nOrigin: %s\nTrigger date: %s\nBarycenter %s",
        event,
        origin,
        date,
        coord2str(center),
    )

    # Connect to the alert database
    matched_alerts: list[Alert] = []
    with dbconnect() as session:
        # Be cautious with yield_per since it may does not allow some querying
        # functionality to work
        for alert in session.query(Alert).order_by(Alert.trigger_date).yield_per(100):
            try:
                alert_loc: Localization = pickle.loads(alert.localization)
                if alert_loc is None:
                    log.debug(
                        "\nEvent: %s\nOrigin: %s\nTrigger date: %s\n NO LOCALIZATION",
                        alert.event,
                        alert.origin,
                        alert.trigger_date,
                    )
                    continue

                # For LVC events barycenter is not that useful info actually, but anyway
                moc_alert = alert_loc.moc()
                alert_center = moc_alert.barycenter()
                log.debug("\nmatching with:")
                log.debug(
                    "\nEvent: %s\nOrigin: %s\nTrigger date: %s\nBarycenter: %s\n",
                    alert.event,
                    alert.origin,
                    alert.trigger_date,
                    coord2str(alert_center),
                )

                # If no time specified, it is considered not a matching
                if alert.trigger_date is None:
                    continue

                # Matched by coordinates
                is_matched_by_coord = not moc_alert.intersection(moc).empty()
                if is_matched_by_coord:
                    log.debug("Events are matched by coordinates")

                # Matched by date and time
                time_diff = abs(Time(date).jd - Time(alert.trigger_date).jd)
                log.debug("Events are located within %.3f days", time_diff)

                is_matched_by_datetime = time_diff <= max_date_diff.to(u.day).value
                if is_matched_by_datetime:
                    log.debug("Events are match by datetime")

                # Are both alerts matching?
                if is_matched_by_coord and is_matched_by_datetime:
                    log.debug("Events are matched")
                    matched_alerts.append(
                        {
                            "id": alert.id,
                            "event": alert.event,
                            "origin": alert.origin,
                            "ra_center": alert.ra_center,
                            "dec_center": alert.dec_center,
                            "error_radius": alert.error_radius,
                            "trigger_date": alert.trigger_date,
                            "localization": alert.localization,
                        }
                    )
            except Exception as e:
                # We do not want the scenario where this function fails on a certain
                # event, and all cross-matching failes with it, so handle all possible
                # exceptions.
                log.error("alerts not crossmatched due to error", exc_info=e)

    if matched_alerts:
        log.debug("Found %d matching alerts.", len(matched_alerts))

    return matched_alerts


def replace_with_matched(
    to_replace: Sequence[Alert], alert_msg: str | bytes, replacement: TargetInfo
):
    log.debug(
        "replacing %d duplicated alerts after cross-matching ...", len(to_replace)
    )
    engine, session = create_session()
    ids = set(a.id for a in to_replace)
    with session:
        session.query(Alert).filter(Alert.id.in_(ids)).delete()
        replacement_alert = Alert()
        replacement_alert.alert_message = (
            bytes(alert_msg, encoding="utf-8")
            if isinstance(alert_msg, str)
            else alert_msg
        )
        replacement_alert.ra_center = replacement.localization.center().ra.deg
        replacement_alert.dec_center = replacement.localization.center().dec.deg
        replacement_alert.error_radius = (
            replacement.localization.error_radius().to_value(u.deg)
        )
        replacement_alert.event = replacement.event
        replacement_alert.origin = replacement.origin
        replacement_alert.importance = replacement.importance
        replacement_alert.trigger_date = replacement.trigger_date
        replacement_alert.localization = pickle.dumps(replacement.localization)
        session.add(replacement_alert)

        try:
            session.commit()
        except Exception as e:
            session.rollback()
            log.error("events not replaced", exc_info=e)


def crossmatch_alerts_by_name(origin: str, event: str) -> list[dict]:
    log.debug("Cross-matching alert with a database ...")
    log.debug("Event %s %s", origin, event)

    matched_alerts: list[Alert] = []
    with dbconnect() as session:
        alerts = session.query(Alert).yield_per(100)
        for alert in alerts:
            log.debug("\nMatching with:")
            log.debug("Event %s %s", alert.origin, alert.event)

            if alert.origin == origin and alert.event == event:
                matched_alerts.append(
                    {
                        "id": alert.id,
                        "event": alert.event,
                        "origin": alert.origin,
                        "ra_center": alert.ra_center,
                        "dec_center": alert.dec_center,
                        "error_radius": alert.error_radius,
                        "trigger_date": alert.trigger_date,
                        "localization": alert.localization,
                    }
                )
                log.debug("Alerts are matching!")

    if matched_alerts:
        log.debug("Found %d matching alerts.", len(matched_alerts))

    return matched_alerts
