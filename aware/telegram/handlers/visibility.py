"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
visibility.py (c) 2024
Desc: Visibility plot of a target
Created:  2024-03-24
Modified: 2024-03-24
"""

from io import BytesIO
import aiogram
from astroplan import FixedTarget
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time

from aware import site
from aware.cache import read_cache, write_cache
from aware.logger import log
from aiogram.utils import markdown
from aiogram.utils.emoji import emojize


async def visibility_handler(message: aiogram.types.Message):
    log.debug(
        "telegram.status(): user asked for visibility plot",
    )
    try:
        cmd_name, name, ra, dec, date, telescope_id = message.text.split(" ")
    except (TypeError, ValueError):
        error_msg = "Invalid message format!"
        return await message.bot.reply_error(error_msg, message)

    try:
        ra = float(ra)
        ra = Angle(ra, unit="deg")
    except ValueError:
        ra = Angle(ra, unit="hr")

    try:
        dec = float(dec)
    except ValueError:
        ...

    dec = Angle(dec, unit="deg")

    if telescope_id not in site.Telescopes:
        await message.bot.reply_error(
            f"Telescope with id={telescope_id} is not available",
            message
        )
    else:
        telescope = site.Telescopes[telescope_id]
        coord = SkyCoord(ra, dec)
        target = FixedTarget(coord, name)
        init_time = Time(date, format="isot")
        key = (
            f"telescope_id={telescope_id} target={name} ra={coord.ra.deg:.0f} "
            f"dec={coord.dec.deg:.0f} date={init_time.datetime.date()}"
        )
        plot_file_content = read_cache(key)
        if plot_file_content is None:
            start_time, end_time = telescope.nearest_observation_window(init_time)
            ax = telescope.plot_airmass(
                [target], start_time=start_time, end_time=end_time
            )
            plot_file = BytesIO()
            ax.get_figure().savefig(plot_file)
            plot_file.seek(0)
            plot_file_content = plot_file.read()
            write_cache(key, plot_file_content)

        plot_file = BytesIO(plot_file_content)

        await message.bot.send_photo(
            message.chat.id,
            aiogram.types.InputFile(plot_file),
            caption="Your visibility plot",
            reply_to_message_id=message.message_id,
        )
