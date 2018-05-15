#!/usr/bin/env python3
"""do some wrap functions"""

import os

from cfutils.align import align
from cfutils.parser import parse_abi, parse_fasta
from cfutils.show import highlight_base, plot_chromatograph


def do_mutation_calling(query_ab1_file, subject_fasta_file):
    """Test plot mutation region"""
    query_record = parse_abi(query_ab1_file, trim=True)
    subject_record = parse_fasta(subject_fasta_file)

    mutations = align(query_record, subject_record, ignore_ambig=True)

    selected_mutation = mutations[3][2]
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(15, 6))
    plot_chromatograph(
        query_record,
        ax,
        xlim=[selected_mutation - 10, selected_mutation + 10])
    highlight_base(selected_mutation, query_record, ax)
    os.makedirs('./temp', exist_ok=True)
    plt.savefig('./temp/test.pdf')
