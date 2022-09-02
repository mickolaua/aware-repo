"""
# ----------------------------------------------------------------------------- 
# Project:     AWARE
# Name:        aware.consumer.py 
# Purpose:     Alert message receiving and handling
# 
# Author:      npank 
# 
# Created:     2022-08-29 
# Copyright:   (c) 2004-2022 AWARE Developers
# ----------------------------------------------------------------------------
"""
from __future__ import annotations

import datetime
from hashlib import md5
from pickle import dumps
from typing import Optional, Sequence

import gcn_kafka
from confluent_kafka import (Consumer, KafkaError, KafkaException, Message,
                             TopicPartition)

from aware.topic import TOPICS

from . import alert, sql
from .logger import log

running = True

# TODO: Create async loop
# FIXME: When timeout is limited, DB continously tries to add the same 
# target after several unique ones.
def consume_loop(
    consumer: gcn_kafka.Consumer,
    topics: Optional[Sequence[str]] = None,
    num_messages: int = 1,
    timeout: float = -1,
    start_date: Optional[datetime.datetime] = None
) -> None:
    global running
    try:
        # No topics provided, assume to subscribe to all available topics
        # cluster_meta = consumer.list_topics()
        # topics = list(cluster_meta.topics.keys())
        log.info("\n" + "-"*79)
        if topics is None:
            topics = TOPICS

        # Calculate offsets from date
        if start_date is not None:
            start_timestamp = int(start_date.timestamp() * 1000)
            partitions = [
                TopicPartition(topic, 0, start_timestamp) for topic in topics]
            offsets = consumer.offsets_for_times(partitions)

            # Assign offset for all topic partitions
            consumer.assign(offsets)

            log.info("started listening for alert messages from %s", start_date)
            
        # consumer.assign(start)
        consumer.subscribe(topics)
        
        # Create SQLite session
        engine, session = sql.create_session()

        # Create alert table if it is not exist (e.g. fresh database)
        if not engine.has_table("alert"):
            sql.Alert.metadata.create_all(engine)
        
        # Consume messages from topics
        while running:
            for message in consumer.consume(num_messages=num_messages, 
                                            timeout=timeout):
                log.info("\n" + " ALERT MESSAGE ".center(79, "#"))
                log.info("\nReceived alert message at %s, checking for errors",
                         datetime.datetime.now())

                # Caught error
                if message.error():
                    err_code = message.error().code()
                    if err_code == KafkaError.UNKNOWN_TOPIC_OR_PART:
                        log.error(
                            "\nTopic %s is not recognized! Go to next message", 
                            message.topic())
                        continue
                else:
                    # Display message in debug mode
                    log.info("No errors found, processing the message")
                    msg: str = message.value()
                    log.debug(msg)

                    # Do not know how to parse alert
                    topic: str = message.topic()
                    if topic not in alert.alert_parsers.parsers:
                        log.warning("Do not know how to parse the alert "
                                    +"message from topic '%s'. Go to next "
                                    +"message", topic)
                        continue

                    # Parse the alert
                    alert_msg = bytes(msg)
                    parser = alert.alert_parsers.get_parser(topic)
                    target_info = parser.parse_alert(alert_msg)
                    log.info("alert was sent by %s on %s as the result of "
                             +" triggering to the event ID %s at %.6f, %+.6f "
                             +"+/- (%.6f, %.6f) deg (stat. sign. %.2f)", 
                             target_info.origin, 
                             target_info.trigger_date, 
                             target_info.event, 
                             target_info.ra_center, 
                             target_info.dec_center, 
                             target_info.error_radius1, 
                             target_info.error_radius2, 
                             target_info.importance)

                    # Unique hash for the alert message to not add duplicates
                    # to the database
                    hash_md5 = md5(alert_msg).hexdigest()

                    # No errors, add the alert to the database
                    localization = dumps(dict(r1=target_info.error_radius1,
                                              r2=target_info.error_radius2))
                    alert_tab = sql.Alert(
                        alert_message=alert_msg, 
                        ra_center=target_info.ra_center, 
                        dec_center=target_info.dec_center,
                        error_radius1=target_info.error_radius1,
                        error_radius2=target_info.error_radius2,
                        localization=localization, 
                        trigger_date=target_info.trigger_date,
                        event=target_info.event, 
                        origin=target_info.origin,
                        importance=target_info.importance,
                        md5=hash_md5)
                    
                    with session:
                        try:
                            session.add(alert_tab)
                        except Exception as e:
                            continue
                        
                        try:
                            session.commit()
                        except Exception as e:
                            log.warning("alert already in the database")
                        else:
                            log.info("alert added to the datebase")
    except KeyboardInterrupt:
        shutdown()

    finally:
        # Close down the consumer to commit final offsets.
        consumer.close()


def shutdown():
    global running
    running = False
