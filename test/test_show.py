#!/usr/bin/env python3

"""test for cfutils"""

import os
import sys
import unittest

from cfutils.parser import parse_abi
from cfutils.show import plot_chromatograph, highlight_base

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError('cfutils requires Python 3.6+!')


class TestShowFunc(unittest.TestCase):
    """Test cfutils/show.py"""

    def test_plot_chromatograph(self):
        """Test plot"""
        from pkg_resources import resource_stream
        input_file = resource_stream(__name__, '../data/A4-3k_SEQ_R_D01.ab1')
        seq = parse_abi(input_file, trim=True)

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        plot_chromatograph(seq, ax, xlim=[22, 44])
        highlight_base(30, seq, ax)
        os.makedirs('./temp', exist_ok=True)
        plt.savefig('./temp/test.pdf')
        input_file.close()


if __name__ == '__main__':
    unittest.main()