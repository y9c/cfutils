#!/usr/bin/env python3
"""test for cfutils."""

import os
import sys
import unittest

from cfutils.parser import parse_abi, parse_fasta
from cfutils.show import highlight_base, plot_chromatograph

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError("cfutils requires Python 3.6+!")


class TestShowFunc(unittest.TestCase):
    """Test cfutils/show.py."""

    def test_plot_chromatograph(self):
        """Test plot."""
        query_record = parse_abi("./data/B5-M13R_B07.ab1")
        subject_record = parse_fasta("./data/3kref.fa")

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        plot_chromatograph(query_record, region=(10, 30), ax=ax)
        highlight_base(14, query_record, ax)
        os.makedirs("./temp", exist_ok=True)
        plt.savefig("./temp/test_plot.pdf")


if __name__ == "__main__":
    unittest.main()
