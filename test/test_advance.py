#!/usr/bin/env python3
"""test for cfutils"""

import os
import sys
import unittest

from cfutils.align import align
from cfutils.parser import parse_abi, parse_fasta
from cfutils.show import highlight_base, plot_chromatograph

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError('cfutils requires Python 3.6+!')


class TestFunc(unittest.TestCase):
    """Test cfutils/align.py"""

    def test_plot_mutation(self):
        """Test plot mutation region"""
        query_record = parse_abi('./data/B5-M13R_B07.ab1')
        subject_record = parse_fasta('./data/3kref.fa')

        mutations = align(query_record, subject_record, ignore_ambig=True)

        selected_mutation = mutations[3][2]
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        plot_chromatograph(
            query_record,
            ax,
            xlim=[selected_mutation - 10, selected_mutation + 10])
        highlight_base(selected_mutation, query_record, ax)
        print(selected_mutation)
        os.makedirs('./temp', exist_ok=True)
        plt.savefig('./temp/test.pdf')


if __name__ == '__main__':
    unittest.main()
