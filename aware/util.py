from __future__ import annotations
import os.path
from uuid import uuid4
from .logger import log
import requests


def download_file(
    url: str,
    timeout: float = 1.0,
    path: str = ".",
    chunk_size: int = 2048,
    overwrite: bool = False,
    raise_file_exist: bool = False
) -> bool:
    status = False

    try:
        filename = os.path.split(url)[1]
    except IndexError as e:
        filename = uuid4().hex
        log.error(
            "could not extract filename from URL: `%s`, assigning unique filename `%s`",
            url,
            filename,
            exc_info=e,
        )

    filename = os.path.join(path, filename)
    if os.path.exists(filename) and os.path.isfile(filename) and not overwrite:
        if raise_file_exist:
            raise FileExistsError(filename)
        else:
            log.warn("file already exists: `%s`, so do not download it", filename)
            status = True
    else:
        try:
            resp = requests.get(url, timeout=timeout, stream=True)
            with open(filename, "wb") as f:
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
            status = True
        except requests.HTTPError as e:
            log.error("could not download file from URL: `%s`", url, exc_info=e)

    return status


def dig_or_sgn_to_super_char_utf8(s: str) -> str:
    """Converts a digit (0-9) or sign char ('-' or '+') to superscript char."""
    if s == "0":
        return "\u2070"
    elif s == "1":
        return "\u00B9"
    elif s == "2":
        return "\u00B2"
    elif s == "3":
        return "\u00B3"
    elif s == "4":
        return "\u2074"
    elif s == "5":
        return "\u2075"
    elif s == "6":
        return "\u2076"
    elif s == "7":
        return "\u2077"
    elif s == "8":
        return "\u2078"
    elif s == "9":
        return "\u2079"
    elif s == "+":
        return "\u207A"
    elif s == "-":
        return "\u207B"
    else:
        raise ValueError(f"character {s} is not a digit or a sign")


def render_number_as_rich_utf8(number: float, precision: int = 6) -> str:
    # For simplicity in further processing just convert the number to a string
    # Right trailing zeros of float numbers are automatically stripped in Python.
    # However, if precision is not specified we want to keep them.
    num_str = f"{number:e}"

    # Split the base and the exponent
    base, exp = num_str.split("e")

    # Clean the exponent
    if set(exp[1:]) == {"0"}:
        exp = ""

    if exp:
        exp_num_part = exp[1:].lstrip("0")
        exp = ("" if exp[0] == "+" else "-") + exp_num_part

    # Split integer and fractional parts of the base
    int_part, frac_part = base.split(".")

    if precision:
        frac_part = str(round(float(f"1.{frac_part}"), precision)).split(".")[1]
        frac_part = frac_part.ljust(precision, "0")
    else:
        frac_part = ""

    # If exponent is just zero return only base
    if exp and precision:
        # if precision:
        new_base = int_part + f".{frac_part}"
        # else:
        #     new_base = int_part
        res = (
            new_base
            + " x 10"
            + "".join(
                [dig_or_sgn_to_super_char_utf8(d) for d in exp]
            )
        )
    else:
        if precision:
            res = int_part + f".{frac_part}"
        else:
            res = f"{int(number):d}"

    return res