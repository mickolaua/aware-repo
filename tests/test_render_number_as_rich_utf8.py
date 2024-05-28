"""
Author: Nicolai Pankov (colinsergesen@gmail.com)
test_render_number_as_rich_utf8 (c) 2024
Desc: test that the number is rendered correctly as rich UTF-8 text
Created:  2024-03-26
Modified: 2024-03-30
"""
from aware.util import render_number_as_rich_utf8


def test_num_without_exponent():
    assert render_number_as_rich_utf8(1.0, 1) == "1.0"
    assert render_number_as_rich_utf8(0.1, 1) == "1.0 x 10\u207B\u00B9"
    assert render_number_as_rich_utf8(1.5, 1) == "1.5"


def test_int_num():
    assert render_number_as_rich_utf8(1, 0) == "1"
    assert render_number_as_rich_utf8(1, 1) == "1.0"
    

def test_num_with_empty_exponent():
    assert render_number_as_rich_utf8(1.00e00, 2) == "1.00"
    assert render_number_as_rich_utf8(0.00e00, 2) == "0.00"


def test_num_precision():
    assert render_number_as_rich_utf8(1503.15, 15) == "1.503150000000000 x 10\u00B3"
    assert render_number_as_rich_utf8(1503, 15) == "1.503000000000000 x 10\u00B3"
    assert render_number_as_rich_utf8(1503, 0) == "1503"
    assert render_number_as_rich_utf8(0.1, 0) == "0"
    assert render_number_as_rich_utf8(1, 0) == "1"
    assert render_number_as_rich_utf8(1.0e-1, 0) == "0"
    assert render_number_as_rich_utf8(0.003, 3) == "3.000 x 10\u207B\u00B3"
    


if __name__ == "__main__":
    test_num_without_exponent()
    test_int_num()
    test_num_with_empty_exponent()


