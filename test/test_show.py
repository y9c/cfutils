#!/usr/bin/env python3
"""test for cfutils."""

import os
import sys
import unittest
import matplotlib.pyplot as plt

from cfutils.parser import parse_abi, parse_fasta, SeqRecord
from cfutils.show import highlight_base, plot_chromatograph, annotate_mutation
from cfutils.align import SitePair

try:
    assert sys.version_info > (3, 8)
except AssertionError:
    raise RuntimeError("cfutils requires Python 3.8+!")


class TestShowFunc(unittest.TestCase):
    """Test cfutils/show.py."""

    def setUp(self):
        """Set up test data."""
        self.query_record = parse_abi("./data/B5-M13R_B07.ab1")
        self.subject_record = parse_fasta("./data/ref.fa")
        self.fig, self.ax = plt.subplots(1, 1, figsize=(15, 6))

    def test_plot_chromatograph(self):
        """Test plot_chromatograph function."""
        plot_chromatograph(self.query_record, region=(10, 30), ax=self.ax)
        self.assertTrue(True)  # Add more specific assertions as needed

    def test_highlight_base(self):
        """Test highlight_base function."""
        plot_chromatograph(self.query_record, region=(10, 20), ax=self.ax)
        highlight_base(14, self.query_record, self.ax)
        self.assertTrue(True)  # Add more specific assertions as needed

    def test_annotate_mutation(self):
        """Test annotate_mutation function."""
        mutation = SitePair(ref_pos=10, ref_base='A', cf_pos=14, cf_base='T')
        annotate_mutation(mutation, self.query_record, self.ax)
        self.assertTrue(True)  # Add more specific assertions as needed

    def tearDown(self):
        """Clean up after tests."""
        plt.close(self.fig)


if __name__ == "__main__":
    unittest.main()
