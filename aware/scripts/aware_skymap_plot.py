import os
from textwrap import dedent
import healpy as hp
from ligo.skymap.plot import cylon
from ligo.skymap.postprocess.util import find_greedy_credible_levels
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.moc import rasterize
import click
from matplotlib import pyplot as plt
import numpy as np
import astropy_healpix as ah
from aware.io import read_hpx_from_moc


@click.command(
    help=dedent(
        """
        Plot skymap and draw the probablity contour on it.
        """
    )
)
@click.argument(
    "input_",
    type=click.Path(file_okay=True, dir_okay=False, exists=True, resolve_path=True),
)
@click.option("-p", "--prob", default=90.0, show_default=True)
@click.option(
    "-n",
    "--name",
    help="Event name",
)
@click.option("-o", "--output", default="", help="Output plot filename")
def main(input_: str, output: str, prob: float, name: str):
    try:
        hpx, hdr = read_hpx_from_moc(input_, header=True)
    except KeyError:
        # Not a MOC skymap
        hpx, hdr = hp.read_map(input_, h=True, nest=True)

    nside = hp.get_nside(hpx)
    c = find_greedy_credible_levels(hpx)

    # Contour area (degrees squared)
    nside = hp.get_nside(hpx)
    pix_area = hp.nside2pixarea(nside, degrees=True)
    area = len(c[c < prob / 100]) * pix_area

    fig = plt.figure(figsize=(10, 6), dpi=100)
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection="astro hours mollweide")
    cs = ax.contour_hpx(
        (100 * c, "ICRS"),
        colors="k",
        linewidths=2,
        levels=[prob],
    )
    plt.clabel(cs, fontsize=6, inline=True)
    title = dedent(
        f"""
        {name if name else "Unknown event"}
        {prob:.1f}% = {area:.0f} deg$^2$
        """
    )
    ax.set_title(title)

    ax.grid(True)
    ax.imshow_hpx((hpx, "icrs"), cmap="cylon")

    if not output:
        output = os.path.splitext(input_)[0] + ".png"
    fig.savefig(output)


if __name__ == "__main__":
    main()