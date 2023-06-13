from matplotlib import pyplot as plt
from aware.visualization.find_chart import plot_find_chart
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from numpy.random import RandomState, random
import os


def test(tmpdir):
    f = tmpdir.mkdir(".tmp").join("chart.png")
    state = RandomState(seed=999)
    ra0 = random(size=1)*360
    dec0 = (random(size=1)*180) - 90
    center = SkyCoord(ra0*u.deg, dec0*u.deg)
    radius = 15*u.arcmin
    ax = plot_find_chart(center, radius, "Random Object")
    ax.figure.savefig(f.basename)

    print(f.basename)

    assert os.path.exists(f.basename), "Find chart not saved!"
    assert os.path.getsize(f.basename) > 0, "Find chart is empty!"


if __name__ == "__main__":
    test()
    