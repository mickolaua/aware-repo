"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
cosmology (c) 2023
Desc: description
Created:  2023-03-05
Modified: 2023-03-05
"""
from __future__ import annotations


from astropy.cosmology import FlatLambdaCDM
from astropy import units as u


# Cosmological constants
H0 = 69.6*u.km/u.s/u.Mpc
Om0 = 0.286

# Flat FRW cosmology
cosmos = FlatLambdaCDM(H0=H0, Om0=Om0)
