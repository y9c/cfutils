#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Chromatogram File Utils.

do some wrap functions
"""

import os
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt

from cfutils.align import align
from cfutils.parser import parse_abi, parse_fasta
from cfutils.show import annotate_mutation, highlight_base, plot_chromatograph


def do_mutation_calling(
    query_ab1_file,
    subject_fasta_file,
    output_dir=None,
    file_basename=None,
    report_mut_info=True,
    report_mut_plot=False,
):
    """test plot mutation region."""
    min_qual = 50
    if not output_dir:
        output_dir = os.path.join(
            os.getcwd(),
            "CFresult_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
        )
    os.makedirs(output_dir, exist_ok=True)

    if not file_basename:
        file_basename = (
            Path(query_ab1_file).stem + "_vs_" + Path(subject_fasta_file).stem
        )

    query_record = parse_abi(query_ab1_file)
    subject_record = parse_fasta(subject_fasta_file)
    mutations = align(query_record, subject_record, ignore_ambig=True)
    # save tsv file
    with open(os.path.join(output_dir, file_basename + ".tsv"), "w") as f_mut:
        f_mut.write(
            "\t".join(
                ["RefLocation", "RefBase", "CfLocation", "CfBase", "CfQual"]
            )
            + "\n"
        )
        for mut in mutations:
            f_mut.write(
                f"{mut.ref_pos}\t{mut.ref_base}\t{mut.cf_pos}\t{mut.cf_base}\t{mut.cf_qual}\n"
            )

    mutations = [m for m in mutations if m.cf_qual and m.cf_qual >= min_qual]
    if mutations and report_mut_plot:
        fig, axes = plt.subplots(
            len(mutations), figsize=(15, 5 * len(mutations))
        )
        flanking_size = 10
        for i, mutation_info in enumerate(mutations):
            if len(mutations) == 1:
                ax = axes
            else:
                ax = axes[i]
            plot_chromatograph(
                query_record,
                ax,
                region=(
                    mutation_info.cf_pos - flanking_size,
                    mutation_info.cf_pos + flanking_size,
                ),
            )
            highlight_base(mutation_info.cf_pos, query_record, ax)
            annotate_mutation(
                [
                    mutation_info.ref_pos,
                    mutation_info.ref_base,
                    mutation_info.cf_pos,
                    mutation_info.cf_base,
                ],
                query_record,
                ax,
            )
        fig.savefig(os.path.join(output_dir, file_basename + ".pdf"))
