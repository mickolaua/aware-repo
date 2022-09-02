**"AWARE"** is acronim for **"Alert Watcher and Astronomical Rapid Exploration"** is a Python 
application for receiving and processing the Kafka Messages distributed via the 
`GCN/TAN <gcn.nasa.gov>`_ broker for alerts on high-energy transient events such as GRBs and 
gravitational wave events detected by LIGO/Virgo/KAGRA interferometers. 
Besides that, **AWARE** provides optimal scheduling for the observations of these astronomical phenomena with telescopes.
The AWARE uses `gcn_kafka <https://github.com/nasa-gcn/gcn-kafka-python>`_ and `confluent_kafka <https://github.com/confluentinc/confluent-kafka-python>`_ under the hood, which provide 
convinient way for working with GCN/TAN Confluent Kafka broker. 

Installation
============
Installation process is quite simple via pip:

``pip install aware``

To install from source:

``pip install .``
