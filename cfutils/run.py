#!/usr/bin/env python3
"""do some wrap functions"""

import os

import matplotlib.pyplot as plt

from cfutils.align import align
from cfutils.parser import parse_abi, parse_fasta
from cfutils.show import plot_chromatograph, highlight_base, annotate_mutation


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
            ["RefLocation", "RefBase", "CfLocation", "CfBase", "CfQual"]) +
                    "\n")
        for m in mutations:
            f_mut.write(
                f"{m.ref_position}\t{m.ref_base}\t{m.cf_position}\t{m.cf_base}\t{m.cf_qual}\n"
            )

    mutations = [m for m in mutations if m.cf_qual >= 50]
    if report_mut_plot:
        fig, ax = plt.subplots(
            len(mutations), figsize=(15, 5 * len(mutations)))
        flanking_size = 10
        for i, mutation_info in enumerate(mutations):
            plot_chromatograph(
                query_record,
                ax[i],
                region=(mutation_info.cf_position - flanking_size,
                        mutation_info.cf_position + flanking_size))
            highlight_base(mutation_info.cf_position, query_record, ax[i])
            annotate_mutation([
                mutation_info.ref_position, mutation_info.ref_base,
                mutation_info.cf_position, mutation_info.cf_base
            ], query_record, ax[i])
        fig.savefig(mut_plot_file)
