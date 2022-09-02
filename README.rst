**"AWARE"** is acronim for **"Alert Watcher and Astronomical Rapid Exploration"** is a Python 
application for receiving and processing the Kafka Messages distributed via the 
`GCN/TAN <gcn.nasa.gov>`_ broker for alerts on high-energy transient events such as GRBs and 
gravitational wave events detected by LIGO/Virgo/KAGRA interferometers. 
Besides that, **AWARE** provides optimal scheduling for the observations of these astronomical phenomena with telescopes.
The AWARE uses `gcn_kafka <https://github.com/nasa-gcn/gcn-kafka-python>`_ and `confluent_kafka <https://github.com/confluentinc/confluent-kafka-python>`_ under the hood, which provide 
convinient way for working with GCN/TAN Confluent Kafka broker. 

Installation
============
Download an .whl file from targs and run pip:

``pip install AWARE-0.1.0-py3-none-any.whl``

Notes on target sorting
=======================

In the current version the automatic target sorting is implemented via script
``aware_sort_targets`` and based on nearest neighbor algorithm. First off, 
algorithm sorts all objects by airmass: inside observational time window, 
the object with highest airmass to time ratio is choosen. Then, for this and 
for all subsequent objects its nearest neighbor is found -- next observational 
target.
