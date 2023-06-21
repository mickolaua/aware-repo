# -*- coding: utf-8 -*-
# @Author: npank
# @Date:   2022-07-21
# @Last Modified by:   npank
# @Last Modified time: 2022-08-29
from aware.credentials import get_credentials
from aware.consumer import prepare_consumer
import pytest
import os


@pytest.mark.skipif(
    not os.getenv("GCN_KAFKA_CLIENT_ID") or not os.getenv("GCN_KAFKA_CLIENT_SECRET"),
    reason="No GCN credentials provided",
)
def test():
    credits = get_credentials()
    consumer = prepare_consumer({}, credits)
    topics = consumer.list_topics().topics
    assert topics is not None, "consumer must return not empty list of topics"


if __name__ == "__main__":
    test()
