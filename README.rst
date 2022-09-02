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

``pip install AWARE-0.1.0-py3-none-any.whl --no-build-isolation``

Receiving GCN alerts
====================
Before receiving alerts, one needs to authorize to `<gcn.nasa.gov>`_ and get credentials (id and secret).
After it's done, set-up enviromental variables (see example for PowerShell):

``$Env:GCN_KAFKA_CLIENT_ID="<your actual id>"``
``$Env:GCN_KAFKA_CLIENT_SECRET="<your actual secret>"``

To run AWARE the YAML config file is needed, default one can be found in this repository. 
Access to the config file is provided via enviromental variable (see example for PowerShell):

``$Env:AWARE_CONFIG_FILE="path/to/config.yaml"``

Now, everything is ready to receive alerts AWARE must be ran from command line:

``python -m aware``

By default, all alerts will be saved to "alert.db" database in the directory from which python is executed.


Target sorting
==============

Optimal observation of targets is needed for fast targeting search of optical transients 
around Glade+ galaxies inside LVC localizaions in the sky. To find the optimal observation 
order of targets, several steps must be followed:

1. Given the desired observational time (probably from alerts on transients),
find the nearest time window, when the site can observe targets (two times: start and end of such period)

2. Select only targets that can be observed in the time window by the site (observer)

3. Sort targets in such order to provided best efficiency of their observation, i.e. transition between targets must be as fast as possible. At the same time, it is desirable to observe as many targets as possible, when transitioning between targets.

Regarding step #3, it is difficult to implement such algorithm. However, the straight forward approach can be used: Nearest Neighbor algorithm.  
The Nearest Neighbor algorithm consists of few steps:

1. Given the time window, find the target that more than other located at airmass > 3 (i.e. airmass-to-time ratio). 
This target will be observed first

2. For this and for all subsequent objects:

  a. find its nearest neighbor
  b. place current target after previous
  c. place nearest neighbor after current target

In the current version, the target sorting is implemented via script ``aware_sort_targets``. 

Usage:

``python aware_sort_targets.py [-h] [-i INPUT] [-o OUTPUT] [-t TIME] [-s SITE] [--airmass-plot AIRMASS_PLOT]``

where,

- ``INPUT`` is the input file storing the list of targets (target name, ra, dec) in SQLITE DB, JSON or CSV format
- ``OUTPUT`` is the name of the file storing the list of sorted targets in JSON format (site coordinates and name, target name, ra, dec)
- ``TIME`` is the desired time when targets should be observed (used as guess for nearest observational time window)
- ``SITE`` short name of the site (observer)
- ``AIRMASS_PLOT`` name of the file with airmass plot (file is not saved if this argument not provided)

