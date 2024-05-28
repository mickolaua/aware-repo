import click
from aware.site import Telescopes
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroplan import FixedTarget
from datetime import datetime
from astropy.time import Time


@click.command()
@click.option(
    "-a",
    "--ra",
    required=True,
    help="Right ascension of the target in astropy notation",
)
@click.option(
    "-d", "--dec", required=True, help="Declination of the target in astropy notation"
)
@click.option(
    "-t",
    "--telescope",
    required=True,
    help="Telescope for which the target visibility should be calculated",
    type=click.Choice(list(Telescopes.keys())),
)
@click.option("-n", "--name", default="", help="name of target")
@click.option(
    "-o", "--output", default="visibility.png", help="output visibility plot file"
)
def main(ra: float, dec: float, telescope: str, name: str, output: str):
    target = FixedTarget(SkyCoord(ra, dec), name=name)
    observer = Telescopes.get(telescope)
    now = Time(datetime.now())
    start_time, end_time = observer.nearest_observation_window(now)
    ax = observer.plot_airmass([target], start_time=start_time, end_time=end_time)
    fig = ax.get_figure()
    fig.savefig(output)
    

if __name__ == "__main__":
    main()
