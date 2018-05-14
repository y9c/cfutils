#!/usr/bin/env python3

"""test for cfutils"""

import sys
import unittest

from cfutils.parser import parse_abi
from cfutils.align import run_align

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError('cfutils requires Python 3.6+!')


class TestAlignFunc(unittest.TestCase):
    """Test cfutils/align.py"""

    def test_plot_chromatograph(self):
        """Test align"""
        from pkg_resources import resource_stream
        input_file = resource_stream(__name__, '../data/A2-3k_SEQ_R_B01.ab1')
        seq = parse_abi(input_file, trim=True)

        subject_fasta = './data/3kref.fa'
        run_align(seq, subject_fasta, ignore_ambiguous=True)
        input_file.close()


if __name__ == '__main__':
    unittest.main()
