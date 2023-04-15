import os
from aware.alert import alert
from io import BytesIO
from lxml import etree
import re
import aiohttp
import asyncio
from astropy.io import fits
import astropy_healpix as ah
import astropy.units as u
from mocpy import MOC
import numpy as np
import pendulum


async def test():
    EVENT_FILE = "tests/gw_event.xml"
    root = alert.get_xml_root(open(EVENT_FILE, "rb").read())

    obs_id = root.find(
        ".//WhereWhen/ObsDataLocation/ObservatoryLocation").attrib["id"]
    if not re.findall(r"(LIGO)|(Virgo)|(LVK)|(LVC)", obs_id):
        raise ValueError

    UT = pendulum.parse(root.find(".//ISOTime").text)
    print(UT)

    for p in root.findall(".//What/Param"):
        if p.attrib["name"] == "GraceID":
            event_id = p.attrib["value"]
    
    print(event_id)

    for p in root.findall(".//What/Group/Param"):
        if p.attrib["name"] == "skymap_fits":
            healpix_url = p.attrib["value"]

        if p.attrib["name"] == "HasNS":
            bns_prob = p.attrib["value"]

    print(healpix_url)
    print(bns_prob)

    for child in root.getchildren():
        if child.tag == "What":
            print(child.getchildren())


    HEALPIX_FNAME = healpix_url.split("/")[-1].strip(",1")
    if not os.path.exists(HEALPIX_FNAME):
        async with aiohttp.request("GET", healpix_url) as r:
            healpix = await r.read()

        with open(HEALPIX_FNAME, "wb") as f:
            f.write(healpix)

    data = fits.open(HEALPIX_FNAME, mode="readonly")[1].data
    uniq = data['UNIQ']
    probdensity = data['PROBDENSITY']
    print(probdensity)

    level, ipix = ah.uniq_to_level_ipix(uniq)
    nside = ah.level_to_nside(level)
    lon, lat = ah.healpix_to_lonlat(ipix, nside, order="nested")
    # loc_dict = 



    # msg = StringIO(etree.tostring(root, pretty_print=True, encoding='unicode'))
    # parser = FermiGBMParser(msg)
    # target_info = parser.parse_alert()

    # assert (
    #     target_info.ra_center == 0.0
    #     and target_info.dec_center == 0.0
    #     and target_info.name == "682918137"
    #     and target_info.origin == "Fermi Satellite, GBM Instrument"
    #     and target_info.error_radius == 0.0
    #     and target_info.time == "2022-08-23T03:28:52.39"
    # ), "Alert was parsed incorrectly"


if __name__ == "__main__":
    asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())
    loop = asyncio.get_event_loop()
    loop.run_until_complete(test())
