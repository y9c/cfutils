#!/usr/bin/env python3
"""test for cfutils"""

import os
import sys
import unittest

from cfutils.align import align
from cfutils.parser import parse_abi
from cfutils.show import highlight_base, plot_chromatograph

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError('cfutils requires Python 3.6+!')


class TestFunc(unittest.TestCase):
    """Test cfutils/align.py"""

    def test_plot_mutation(self):
        """Test plot mutation region"""
        from pkg_resources import resource_stream
        input_file = resource_stream(__name__, '../data/B5-M13R_B07.ab1')
        seq = parse_abi(input_file, trim=True)

        subject_fasta = './data/3kref.fa'
        mutations = align(seq, subject_fasta, ignore_ambig=True)
        input_file.close()

        #  selected_mutation = mutations[5][2]
        selected_mutation = 224
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        plot_chromatograph(
            seq, ax, xlim=[selected_mutation - 10, selected_mutation + 10])
        highlight_base(selected_mutation, seq, ax)
        print(selected_mutation)
        os.makedirs('./temp', exist_ok=True)
        plt.savefig('./temp/test.pdf')
        input_file.close()


if __name__ == '__main__':
    unittest.main()
