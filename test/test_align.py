#!/usr/bin/env python3
"""test for cfutils."""

import sys
import unittest

from cfutils.align import align_chromatograph
from cfutils.parser import parse_abi

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError("cfutils requires Python 3.6+!")


class TestAlignFunc(unittest.TestCase):
    """Test cfutils/align.py."""

    def test_align(self):
        """Test align."""
        from pkg_resources import resource_stream

        input_file = resource_stream(__name__, "../data/C5-3k_SEQ_R_A04.ab1")
        seq = parse_abi(input_file)

        subject_fasta = "./data/3kref.fa"
        align_chromatograph(seq, subject_fasta, ignore_ambig=True)
        input_file.close()


if __name__ == "__main__":
    unittest.main()
