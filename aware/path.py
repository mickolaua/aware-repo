"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
path.py (c) 2024
Desc: Miscellaneous utilities for working with path
Created:  2024-05-28
Modified: !date!
"""
import os


def get_full_path(path: str):
    full_path = os.path.realpath(os.path.expanduser(path))
    return full_path
