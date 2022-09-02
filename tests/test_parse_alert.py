from aware import alert
from io import StringIO
from lxml import etree
from aware.logger import log


def test():
    root = etree.parse("New Text Document.xml")

    print(alert.alert_parsers)

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
    test()