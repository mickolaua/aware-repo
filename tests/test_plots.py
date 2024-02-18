from astroplan.target import FixedTarget
from astropy.coordinates import Angle, SkyCoord
from matplotlib import pyplot as plt
import pytest
from aware.site.main import Abastumani_AS32
from astropy.time import Time, TimeDelta
from aware.alert.plugins.lvc import LVC_PRELIMINARY_Parser
from aware.visualization import plot_moc
from aware.visualization.find_chart import plot_find_chart
import os
import matplotlib as mpl
from astropy import units as u


def test_vis_plot(tmpdir: str):
    mpl.use("agg")
    c = SkyCoord(120, 50, unit=["deg"] * 2)
    t = FixedTarget(c, name="A target")
    fn = tmpdir.join("plot.png")
    Abastumani_AS32.plot_airmass(
        [t],
        start_time=Time("2023-03-19T23:00:00", format="isot"),
        end_time=Time("2023-03-20T02:00:00"),
        fig_path=fn,
    )

    assert os.path.exists(fn) and os.path.getsize(
        fn
    ), "Not empty plot file should exist"



def test_find_chart_plot(tmpdir: str):
    mpl.use("agg")
    c = SkyCoord(0, 0, unit=["deg"] * 2)
    R = Angle(15*u.arcmin)
    fn = tmpdir.join("plot.png")
    ax = plot_find_chart(c, R, "")
    ax.get_figure().savefig(fn)
    assert os.path.exists(fn) and os.path.getsize(
        fn
    ), "Not empty plot file should exist"


@pytest.mark.skip(reason="LVC skymap is not available anymore")
def test_lvk_plot(test_dir: str, tmpdir: str):
    mpl.use("agg")
    with open(os.path.join(test_dir, "alert_messages", "gw_event.xml"), "rb") as f:
        msg = f.read()
    info = LVC_PRELIMINARY_Parser.parse_alert(msg)
    moc = info.localization.moc()
    ax = info.localization.plot(
        Abastumani_AS32,
        Time("2023-03-19T23:00:00", format="isot"),
        Time("2023-03-20T02:00:00"),
        [],
    )
    fn = tmpdir.join("plot.png")
    ax.figure.savefig(fn)

    assert os.path.exists(fn) and os.path.getsize(
        fn
    ), "Not empty plot file should exist"


if __name__ == "__main__":
    test_vis_plot(".tmp")
    test_lvk_plot("tests", ".tmp")
