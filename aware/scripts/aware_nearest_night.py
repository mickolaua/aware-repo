from aware.site.main import Telescopes
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroplan.target import FixedTarget
from datetime import datetime
from astropy.time import Time
from matplotlib import pyplot as plt
import click
from aware.logger import log



@click.command()
@click.argument("telescope", type=str)
def cli(telescope: str):
	t = Telescopes[telescope]
	date = Time(datetime.now())
	start, stop = t.nearest_observation_window(date)
	log.info("Nearest night at %s: since %s to %s.", t.full_name,
		start.isot, stop.isot)
	
	print(start.isot)
	print(stop.isot)


def main():
	cli()


if __name__ == '__main__':
	main()
