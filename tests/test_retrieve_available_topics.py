# -*- coding: utf-8 -*-
# @Author: npank
# @Date:   2022-07-21
# @Last Modified by:   npank
# @Last Modified time: 2022-08-29
from gcn_kafka import Consumer
from aware.config import CLIENT_ID, CLIENT_SECRET

def test():
	consumer = Consumer(client_id=CLIENT_ID, client_secret=CLIENT_SECRET)	
	topics = consumer.list_topics().topics

	assert topics is not None, "consumer must return not empty list of topics"

if __name__ == '__main__':
	test()