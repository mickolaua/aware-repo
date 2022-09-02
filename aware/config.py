# -*- coding: utf-8 -*-
# @Author: npank
# @Date:   2022-07-21
# @Last Modified by:   npank
# @Last Modified time: 2022-08-29
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Mapping
import yaml
import dotenv

ABS_DIR = os.path.abspath(os.path.dirname(__file__))

AWARE_CONFIG_FILE = os.getenv("AWARE_CONFIG_FILE", None)
if not AWARE_CONFIG_FILE:
    # Config file is not found
    raise FileNotFoundError("Config file not found.")


def load_config(
    path: str | Path, loader: yaml.Loader = yaml.SafeLoader
) -> Mapping[str, Any]:
    with open(os.path.expanduser(os.path.realpath(path)), "rb") as f:
        cfg = yaml.load(f, yaml.SafeLoader)

    return cfg


# Parse config
cfg = load_config(AWARE_CONFIG_FILE)

if cfg["dev"]:
    if not dotenv.load_dotenv(verbose=True):
        raise FileNotFoundError(".env file is not found")

# Credentials
try:
    CLIENT_ID = os.environ["GCN_KAFKA_CLIENT_ID"]
except KeyError as e:
    raise RuntimeError("client id for GCN Kafka broker is not set")

try:
    CLIENT_SECRET = os.environ["GCN_KAFKA_CLIENT_SECRET"]
except KeyError as e:
    raise RuntimeError(
        "client secret for GCN Kafka broker is not specified"
    )

