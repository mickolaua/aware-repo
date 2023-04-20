from astropy.coordinates import SkyCoord


def coord2str(coord: SkyCoord) -> tuple[str, str]:
    h, m, s = coord.ra.hms
    s_rem = int((s - int(s)) * 10)
    ra = "{:02d}:{:02d}:{:02d}.{:01d}".format(int(h), int(m), int(s), s_rem)
    d, m, s = coord.dec.dms
    dec = "{:+02d}:{:02d}:{:02d}".format(int(d), abs(int(m)), abs(int(s)))

    return ra, dec
