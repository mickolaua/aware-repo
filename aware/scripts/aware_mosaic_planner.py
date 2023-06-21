from textwrap import dedent
import healpy as hp
from healpy.newvisufunc import projview
from healpy.visufunc import cartview, mollview
from functools import partial
import json
import os
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy.ma as ma
import shutil
import click
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.moc import rasterize
import astropy_healpix as ah
from aware.logger import log
from aware.site import Telescopes
from aware.field import Field
from aware.planning.program import create_observation_program
from astropy.coordinates import SkyCoord
from astropy import units as u
from aware.localization.main import LVCSkyMap
from ligo.skymap.plot import cylon
from ligo.skymap.postprocess.util import find_greedy_credible_levels
from aware.io import read_hpx_from_moc
from mocpy import MOC


def plot_hpx(hpx, targets: list[Field]) -> Axes:
    c = find_greedy_credible_levels(hpx)

    fig = plt.figure(figsize=(10, 6), dpi=100)
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], projection="astro hours mollweide")
    cs = ax.contour_hpx(
        (100 * c, "ICRS"),
        colors="k",
        linewidths=2,
        levels=[90],
    )
    plt.clabel(cs, fontsize=6, inline=True)

    ax.grid(True)
    ax.imshow_hpx((hpx, "icrs"), cmap="cylon")

    for f in targets:
        box = f.box()
        rect = Rectangle(
            (box.left_lower.ra.deg, box.left_lower.dec.deg),
            width=f.width.to_value("deg"),
            height=f.height.to_value("deg"),
            transform=ax.get_transform("world"),
            edgecolor="k",
            lw=3,
            alpha=0.3,
        )
        ax.add_artist(rect)

    return ax


def cur_matrix(
    arr, y_n, x_n, scale_x, scale_y, alpha=10 ** (5), beta=10 ** (-2), max_step=5
):
    cur_arr = np.zeros_like(arr)
    bord1 = [int(y_n - max_step) if int(y_n - max_step) > 0 else 0][0]
    bord2 = [int(y_n + max_step) if int(y_n + max_step) < scale_y else int(scale_y)][0]
    for i in range(bord1, bord2):
        for j in range(scale_x):
            dist_y = abs(y_n - i)
            abs_x = abs(x_n - j)
            dist_x = [abs_x if (abs_x <= (scale_x / 2)) else (scale_x - abs_x)][0]
            if dist_x > max_step:
                cur_arr[i][j] = 0
            else:
                cur_arr[i][j] = (alpha * arr[i][j]) - (
                    beta * (dist_y**2 + dist_x**2) ** 0.5
                )
    return cur_arr


def simple(arr, i, scale_x, scale_y):
    report = []
    pred_xy = []
    y_n, x_n = np.unravel_index(np.argmax(arr), arr.shape)
    val_n = arr[y_n][x_n]
    report = [y_n, x_n, f"{int(val_n)}%"]
    arr[y_n][x_n] = 0
    pred_xy.append(report)
    for z in range(i - 1):
        cur_arr = cur_matrix(arr, y_n, x_n, scale_x, scale_y)
        y_n, x_n = np.unravel_index(np.argmax(cur_arr), cur_arr.shape)
        val_n = arr[y_n][x_n]
        report = [y_n, x_n, f"{int(val_n)}%"]
        arr[y_n][x_n] = 0
        pred_xy.append(report)
    return pred_xy


def transition_xy(pred_xy, scale_x, scale_y):
    for i in range(len(pred_xy)):
        cur_y = pred_xy[i][0] + 0.5
        cur_x = pred_xy[i][1] + 0.5
        d_ang_y = 180
        d_ang_x = 360
        pred_xy[i][0] = round(((d_ang_y / scale_y) * cur_y) - 90, 2)
        pred_xy[i][1] = round(
            (((d_ang_x / scale_x) * (scale_x - cur_x)) + 180) % 360, 2
        )
    return pred_xy


def intt(num):
    num = int(num + (0.5 if num > 0 else -0.5))
    return num


def write_json(path, data):
    dirname = os.path.dirname(path)
    os.makedirs(dirname, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=4, sort_keys=True)


def read_json(target_path: str) -> dict:
    with open(target_path, "r") as f:
        data = json.load(f)

    return data


def get_list(
    hpx,
    file_name,
    day,
    matrix_ang,
    observatory_name,
    frame_longitude,
    frame_latitude,
    iteration,
    border_longitude,
    border_latitude,
):
    decay = 0.9
    npix = len(hpx)
    nside = hp.pixelfunc.get_nside(hpx)
    frame_x = decay * frame_longitude
    frame_y = decay * frame_latitude
    scale_x = 360 / frame_x
    scale_y = 180 / frame_y
    scale_x = intt(scale_x)
    scale_y = intt(scale_y)
    border_lat = np.empty_like(border_latitude)
    border_long = np.empty_like(border_longitude)
    d_ang_y = 180
    d_ang_x = 360

    border_lat[0] = (border_latitude[0] + 90) * scale_y / d_ang_y
    border_lat[1] = (border_latitude[1] + 90) * scale_y / d_ang_y
    border_long[0] = [
        -((border_longitude[0] - 180) * scale_x / d_ang_x)
        if border_longitude[0] < 180
        else scale_x - ((border_longitude[0] - 180) * scale_x / d_ang_x)
    ][0]
    border_long[1] = [
        -((border_longitude[1] - 180) * scale_x / d_ang_x)
        if border_longitude[1] < 180
        else scale_x - ((border_longitude[1] - 180) * scale_x / d_ang_x)
    ][0]
    proj = hp.projector.CartesianProj(hpx, xsize=scale_x, ysize=scale_y)
    map = proj.projmap(map=hpx, vec2pix_func=partial(hp.vec2pix, nside))
    arr = map.copy()
    mx = np.amax(arr)
    arr = (arr / mx) * 100

    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if i < border_lat[0] or i > border_lat[1]:
                arr[i][j] = 0
            if border_long[0] > border_long[1]:
                if j > border_long[0] or j < border_long[1]:
                    arr[i][j] = 0
            if border_long[0] < border_long[1]:
                if border_long[0] < j < border_long[1]:
                    arr[i][j] = 0
    pred_xy = []
    pred_ang = []
    pred_xy = simple(arr, iteration, scale_x, scale_y)
    pred_ang = transition_xy(pred_xy, scale_x, scale_y)
    day = f"{day}_day"
    root_dir = os.path.dirname(file_name)
    outdir = os.path.join(root_dir, day, observatory_name)
    write_json(os.path.join(outdir, observatory_name + ".json"), pred_ang)
    scope = Telescopes[observatory_name]
    RA = np.asarray(pred_ang)[:, 1].astype(float)
    DEC = np.asarray(pred_ang)[:, 0].astype(float)
    fields = [
        Field(
            center_coord=SkyCoord(ra * u.deg, dec * u.deg),
            width=scope.fov.height,
            height=scope.fov.height,
        )
        for ra, dec in zip(RA, DEC)
    ]

    program = create_observation_program(scope, fields)
    if program:
        with open(os.path.join(outdir, observatory_name + ".list"), "w") as f:
            f.write(program)

    beta = 1 / decay
    theta_arr = []
    phi_arr = []
    for j in range(len(hpx)):
        theta, phi = matrix_ang[j]
        for i in range(len(pred_ang)):
            my_theta, my_phi, trash = pred_ang[i]
            if (
                my_theta - ((1 / decay) * frame_y / 2)
                < theta
                < my_theta + ((1 / decay) * frame_y / 2)
            ) and (
                my_phi - ((1 / decay) * frame_x / 2)
                < phi
                < my_phi + ((1 / decay) * frame_x / 2)
            ):
                theta_arr.append(my_theta)
                phi_arr.append(my_phi)
                hpx[j] = 0

    # Specify the dtype to shut down the Healpy message on type changed
    hp.write_map(file_name, hpx, overwrite=True, dtype=np.float64)

    draw = True
    if draw == True:
        # projview(
        #     hpx,
        #     coord=["G"],
        #     title=f"{observatory_name} with {iteration} iter and borders long: {border_longitude}, lat: {border_latitude}",
        #     width=10,
        #     graticule=True,
        #     graticule_labels=True,
        #     unit="cbar label",
        #     xlabel="longitude",
        #     ylabel="latitude",
        #     cb_orientation="horizontal",
        #     projection_type="mollweide",
        # )
        # plt.plot(phi_arr, theta_arr, ls="", marker="s", color="r", ms=10)
        # plt.gca().invert_xaxis()
        ax = plot_hpx(hpx, fields)
        dirname = os.path.dirname(file_name)
        ax.get_figure().savefig(
           os.path.join(outdir, f"skymap_{observatory_name}_{iteration}.png")
        )


@click.command(
    help=dedent(
        """
        Perform day-by-day skymap mosaic planning for the telescopes specified in the 
        config.
        """
    )
)
@click.option(
    "-w",
    "--widefield-only",
    is_flag=True,
    default=True,
    help="Plan mosaic for widefield telescopes only",
)
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(file_okay=True, dir_okay=False, exists=True, resolve_path=True),
    help="HEALPix skymap file path",
)
@click.option(
    "-p",
    "--program-file",
    required=True,
    type=click.Path(file_okay=True, dir_okay=False, exists=True, resolve_path=True),
    help="Path to JSON program file",
)
def main(input: str, widefield_only: bool, program_file: str):
    decay = 0.9  # decay frame coeff
    sky_area = int(4 * 180**2 / np.pi) + 1
    root_dir, orig_file_name = os.path.split(input)
    program = read_json(program_file)
    file_name = os.path.splitext(orig_file_name)[0] + "_curr.fits"

    day = program["day"]
    if day == 1:
        shutil.copyfile(
            os.path.join(root_dir, orig_file_name), os.path.join(root_dir, file_name)
        )
    # hpx = hp.read_map(os.path.join(root_dir, file_name))
    # NOTE: Modern fits skymap are in MOC ordering
    hpx = read_hpx_from_moc(os.path.join(root_dir, file_name))
    nside = hp.get_nside(hpx)

    # Only 90% area is needed
    c = find_greedy_credible_levels(hpx)
    hpx_in_c = np.arange(len(hpx))[c < 0.9]
    new_hpx = hpx.copy()
    new_hpx[~(c < 0.9)] = 0.0
    hpx = new_hpx

    npix = len(hpx)
    matrix_ang = []
    for i in range(len(hpx)):
        cur_val = hpx[i]
        theta, phi = hp.pix2ang(nside, i)
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)
        matrix_ang.append([dec, ra])

    for i in range(len(program["annotations"])):
        observatory_name = program["annotations"][i]["observatory_name"]

        # Check that observatory is supported by AWARE
        supported = observatory_name in Telescopes
        if not supported:
            log.warning("Telescope with id `%s` is not supported", observatory_name)
            continue

        scope = Telescopes[observatory_name]

        # Check that scope is widefield
        if widefield_only and not scope.fov.is_widefield:
            log.warning("Telescope with id `%s` is not widefield", observatory_name)
            continue

        frame_longitude = scope.fov.width.to_value("deg")
        frame_latitude = scope.fov.height.to_value("deg")

        # frame_longitude = config["annotations"][i]["frame_longitude"]
        # frame_latitude = config["annotations"][i]["frame_latitude"]
        iteration = program["annotations"][i]["iteration"]
        border_longitude = program["annotations"][i]["border_longitude"]
        border_latitude = program["annotations"][i]["border_latitude"]
        day = day
        get_list(
            hpx,
            os.path.join(root_dir, file_name),
            day,
            matrix_ang,
            observatory_name,
            frame_longitude,
            frame_latitude,
            iteration,
            border_longitude,
            border_latitude,
        )


if __name__ == "__main__":
    main()
