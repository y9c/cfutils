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
from typing import List

import matplotlib.pyplot as plt

from cfutils.align import call_mutations
from cfutils.parser import parse_abi, parse_fasta
from cfutils.show import annotate_mutation, highlight_base, plot_chromatograph


def do_mutation_showing(
    query_record, mutations: List, output_fig_file: str
) -> None:
    """report mutations in pdf format."""
    min_base_qual = 50
    min_local_qual = 20

    mutations = sorted(mutations, key=lambda x: x.cf_pos)
    flanking_size = 6
    windows_size = 30
    mutation_windows = []
    start_pos = max(1, mutations[0].cf_pos - flanking_size)
    mutation_region = []
    for idx, mut in enumerate(mutations):
        if mut.cf_pos + flanking_size <= start_pos + windows_size:
            mutation_region.append(mut)
        else:
            mutation_windows.append(mutation_region)
            start_pos = max(1, mutations[idx].cf_pos - flanking_size)
            mutation_region = [mut]
    mutation_windows.append(mutation_region)

    fig, axes = plt.subplots(
        len(mutation_windows), figsize=(20, 5 * len(mutation_windows))
    )
    for idx, mutation_region in enumerate(mutation_windows):
        if len(mutation_windows) == 1:
            ax = axes
        else:
            ax = axes[idx]
        region_start = max(1, mutation_region[0].cf_pos - flanking_size)
        plot_chromatograph(
            query_record,
            region=(region_start, region_start + windows_size),
            ax=ax,
        )
        for mut in mutation_region:
            base_passed = (
                mut.qual_site is not None
                and mut.qual_site >= min_base_qual
                and mut.qual_local is not None
                and mut.qual_local >= min_local_qual
            )
            highlight_base(
                mut.cf_pos, query_record, ax, passed_filter=base_passed
            )
            annotate_mutation(mut, query_record, ax)
    fig.savefig(output_fig_file, bbox_inches="tight")


def report_mutation(
    query_ab1_file,
    subject_fasta_file,
    output_dir=None,
    file_basename=None,
    report_mut_info=True,
    report_mut_plot=False,
):
    """reprot mutation within region."""
    if output_dir is None:
        output_dir = os.path.join(
            os.getcwd(),
            "CFresult_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
        )
    os.makedirs(output_dir, exist_ok=True)

    if file_basename is None:
        file_basename = (
            Path(query_ab1_file).stem + "_vs_" + Path(subject_fasta_file).stem
        )

    query_record = parse_abi(query_ab1_file)
    subject_record = parse_fasta(subject_fasta_file)
    mutations = call_mutations(query_record, subject_record, ignore_ambig=True)
    # save tsv file
    with open(os.path.join(output_dir, file_basename + ".tsv"), "w") as f_mut:
        header = [
            "RefLocation",
            "RefBase",
            "CfLocation",
            "CfBase",
            "SiteQual",
            "LocalQual",
        ]
        f_mut.write("\t".join(header) + "\n")
        for mut in mutations:
            f_mut.write(
                f"{mut.ref_pos}\t{mut.ref_base}\t{mut.cf_pos}\t{mut.cf_base}\t{mut.qual_site}\t{mut.qual_local}\n"
            )
    if mutations:
        output_fig_file = os.path.join(output_dir, file_basename + ".pdf")
        do_mutation_showing(query_record, mutations, output_fig_file)
