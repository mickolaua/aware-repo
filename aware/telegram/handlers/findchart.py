"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
findchart.py (c) 2024
Desc: description
Created:  2024-03-25
Modified: 2024-03-25
"""
# @dp.message_handler(commands=["findchart"])
from io import BytesIO
import aiogram

from aware.cache import read_cache, write_cache
from aware.logger import log
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u

from aware.visualization.find_chart import plot_find_chart


async def findchart_handler(message: aiogram.types.Message):
    log.debug(
        "user asked for finding chart",
    )
    try:
        _, name, ra, dec = message.text.split(" ")
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

        coord = SkyCoord(ra, dec)
        key = f"ra={coord.ra.deg:.0f}, dec={coord.dec.deg:.0f}, name={name}"
        plot_file_content = read_cache(key)
        if plot_file_content is None:
            ax = plot_find_chart(coord, 15 * u.arcmin, name)
            plot_file = BytesIO()
            ax.get_figure().savefig(plot_file, format="png")
            plot_file.seek(0)
            plot_file_content = plot_file.read()
            write_cache(key, plot_file_content)

        plot_file = BytesIO(plot_file_content)

        await message.bot.send_photo(
            message.chat.id,
            aiogram.types.InputFile(plot_file),
            reply_to_message_id=message.message_id,
        )
    except IndexError as e:
        log.error(
            "field is not available in SkyView or no connection: %s",
            e,
            exc_info=e,
        )
        await message.bot.reply_error("Field is not available in Sky View!", message)
    except Exception as e:
        log.error(
            "can not create finding chart, message was %s",
            message.text,
            exc_info=e,
        )
        await message.bot.reply_error(
            "Unable to create the finding chart with provided parameters!",
            message
        )