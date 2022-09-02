# -*- coding: utf-8 -*-
# @Author: npank
# @Date:   2022-07-22
# @Last Modified by:   npank
# @Last Modified time: 2022-08-29
from aware.logger import log
from aware.config import cfg
import sys
import os


def test() -> int:
    cfg["aware"]["logger"]["filename"] = "test1.log"
    log.critical("critical output (long)" + "ara-"*79 + "ara")
    log.error("error output")
    log.warning("warning output")
    log.info("info output")
    log.debug("debug output")

    assert open("test1.log", "r").read() is not None


if __name__ == '__main__':
    sys.exit(test())