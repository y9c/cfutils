#!/usr/bin/env python3
"""test for cfutils."""

import sys
import unittest

#  from cfutils.align import align
#  from cfutils.show import highlight_base, plot_chromatograph
#  from cfutils.parser import parse_abi, parse_fasta
from cfutils.run import report_mutation

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError("cfutils requires Python 3.6+!")


class TestFunc(unittest.TestCase):
    """Test cfutils/align.py."""

    def test_plot_mutation(self):
        """Test plot mutation region."""
        report_mutation(
            query_ab1_file="./data/B5-M13R_B07.ab1",
            subject_fasta_file="./data/ref.fa",
            output_dir=None,
            file_basename=None,
            report_all_sites=True,
            report_mut_plot=True,
        )


if __name__ == "__main__":
    unittest.main()
