#!/usr/bin/env python3
"""test for cfutils"""

import os
import sys
import unittest

#  from cfutils.align import align
#  from cfutils.show import highlight_base, plot_chromatograph
#  from cfutils.parser import parse_abi, parse_fasta
from cfutils.run import do_mutation_calling

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError('cfutils requires Python 3.6+!')


class TestFunc(unittest.TestCase):
    """Test cfutils/align.py"""

    def test_plot_mutation(self):
        """Test plot mutation region"""
        do_mutation_calling(
            query_ab1_file='./data/B5-M13R_B07.ab1',
            subject_fasta_file='./data/3kref.fa',
            report_mut_info=True,
            report_mut_plot=True,
            mut_info_file="./temp/test.tsv",
            mut_plot_file="./temp/test.pdf")


if __name__ == '__main__':
    unittest.main()
