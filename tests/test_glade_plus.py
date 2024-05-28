import numpy as np
import pytest
from aware.glade import GladeCatalog
from astropy.coordinates import SkyCoord
from astropy import units as u
from aware.alert.plugins.swift import SWIFTXRTAlertParser
import os.path
import mhealpy as mhp
from astropy.table import Table
from ligo.skymap.io.fits import write_sky_map
from astropy.coordinates import SkyCoord
from mocpy import MOC
from ligo.skymap.moc import uniq2ang


@pytest.mark.skip(reason="Long test and requires local GLADE+ catalog.")
def test():
    hpx = mhp.HealpixBase(order=8, scheme="nested")
    ra0 = 100  # deg
    dec0 = 30  # deg
    radius = 0.2  # deg
    disc_pix = hpx.query_disc(mhp.ang2vec(ra0, dec0, lonlat=True), radius)
    m = mhp.HealpixMap.moc_from_pixels(hpx.nside, disc_pix, density=True, nest=True)

    for pix in np.arange(m.npix):
        ra, dec = m.pix2ang(pix, lonlat=True)
        m[pix] = np.exp(-((ra - ra0) ** 2 + (dec - dec0) ** 2) / 2 / (radius) ** 2)

    uniq = m.uniq
    prob = m.data

    moc = Table({"UNIQ": uniq, "PROBDENSITY": prob})

    galaxies = GladeCatalog.query_skymap_local(
        moc, dist_lo_bound=35 * u.Mpc, dist_hi_bound=40 * u.Mpc
    )

    assert len(galaxies) == pytest.approx(3000, abs=1000)


if __name__ == "__main__":
    test()
