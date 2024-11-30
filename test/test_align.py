#!/usr/bin/env python3
"""test for cfutils."""

import sys
import unittest

from cfutils.align import align_chromatograph, call_mutations
from cfutils.parser import parse_abi, parse_fasta

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError("cfutils requires Python 3.6+!")


class TestAlignFunc(unittest.TestCase):
    """Test cfutils/align.py."""

    def setUp(self):
        self.query_record = parse_abi("./data/B5-M13R_B07.ab1")
        self.subject_record = parse_fasta("./data/ref.fa")

    def test_align_chromatograph(self):
        """Test align_chromatograph function."""
        sitepairs = align_chromatograph(self.query_record, self.subject_record)
        self.assertIsInstance(sitepairs, list)

    def test_call_mutations(self):
        """Test call_mutations function."""
        mutations = call_mutations(self.query_record, self.subject_record)
        self.assertIsInstance(mutations, list)


if __name__ == "__main__":
    unittest.main()
