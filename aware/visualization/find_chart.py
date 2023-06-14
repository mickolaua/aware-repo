from __future__ import annotations

import numpy as np
from astroplan.target import FixedTarget
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io.fits import HDUList, ImageHDU
from astropy.visualization import (
    ImageNormalize,
    PercentileInterval,
    SqrtStretch,
    imshow_norm,
)
from astropy.wcs import WCS
from astroquery.skyview import SkyView
from ligo.skymap.plot.marker import reticle
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib_scalebar.scalebar import ScaleBar


def plot_find_chart(center: SkyCoord, radius: Angle, name: str) -> Axes:
    img: HDUList = SkyView.get_images(
        center,
        survey="DSS",
        coordinates="J2000",
        radius=radius,  # type: ignore
    )[
        0
    ]  # type: ignore
    hdr = img[0].header  # type: ignore
    data = img[0].data  # type: ignore
    vmin = data.mean() - data.std()  # type: ignore
    vmax = data.mean() + data.std()  # type: ignore
    wcs = WCS(hdr)
    ax = plt.subplot(projection=wcs)
    interval = PercentileInterval(90)
    stretch = SqrtStretch()
    imshow_norm(data, ax, interval=interval, stretch=stretch)  # type: ignore
    ax.plot(
        center.ra.deg,  # type: ignore
        center.dec.deg,  # type: ignore
        ls="",
        # Left Top crosshair
        marker=reticle(which="lt"),
        ms=5 * radius.to_value("arcmin"),  # type: ignore
        markeredgewidth=3,
        color="red",
        transform=plt.gca().get_transform("world"),
    )
    # print(wcs.wcs.cdelt*3600)
    scale_bar = ScaleBar(
        wcs.wcs.cdelt[0] * 60,
        units="'",
        dimension="angle",
        scale_loc="top",
        location="lower center",
        height_fraction=0.005,
        frameon=False,
        font_properties={"size": 12, "weight": "bold"},
        scale_formatter=lambda value, dimension: f"{value} arcsec",
    )
    scale_bar.set_color("red")
    ax.add_artist(scale_bar)
    ax.grid(True)
    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")
    ax.legend(
        markerscale=0.4,
        fontsize=14,
        frameon=False,
        labelcolor="red",
        edgecolor="b",
    )
    ax.set_title(f"Finding chart of {name}")

    return ax