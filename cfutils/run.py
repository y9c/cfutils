#!/usr/bin/env python3
"""do some wrap functions"""

import os

import matplotlib.pyplot as plt

from cfutils.align import align
from cfutils.parser import parse_abi, parse_fasta
from cfutils.show import highlight_base, plot_chromatograph


def do_mutation_calling(query_ab1_file,
                        subject_fasta_file,
                        report_mut_info=True,
                        report_mut_plot=False,
                        mut_info_file="./temp/test.tsv",
                        mut_plot_file="./temp/test.pdf"):
    """Test plot mutation region"""
    os.makedirs('./temp', exist_ok=True)

    query_record = parse_abi(query_ab1_file)
    subject_record = parse_fasta(subject_fasta_file)
    mutations = align(query_record, subject_record, ignore_ambig=True)
    # save tsv file
    with open(mut_info_file, "w") as f_mut:
        f_mut.write("\t".join(
            ["RefLocation", "RefBase", "CfLocation", "MutBase"]) + "\n")
        for mutation_info in mutations:
            f_mut.write("\t".join(str(x) for x in mutation_info) + "\n")

    # save pdf file
    #  fig, ax = plt.subplots(3, -(-len(mutations) // 3), figsize=(15, 6))
    # don't forget to use ax in form ax[1, 2]
    #  mutations = mutations[:-5]
    if report_mut_plot:
        fig, ax = plt.subplots(len(mutations), figsize=(15, 5 * len(mutations)))
        flanking_size = 10
        # bug mutation location in cf file out of range
        for i, mutation_info in enumerate(mutations):
            plot_chromatograph(
                query_record,
                ax[i],
                xlim=[
                    mutation_info[2] - flanking_size,
                    mutation_info[2] + flanking_size
                ])
            highlight_base(mutation_info[2], query_record, ax[i])
        fig.savefig(mut_plot_file)
