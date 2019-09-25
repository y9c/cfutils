#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Chromatogram File Utils.

show alignment with matplotlib

author:     Fabio Zanini
date:       09/12/13
content:    Plot functions for Sanger chromatographs.
modified:   By Ye Chang in 2018-05-14
"""

from collections import defaultdict
from typing import Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord

from .align import SitePair, align_chromatograph
from .utils import get_logger

LOGGER = get_logger(__name__)


def plot_chromatograph(
    seq: SeqRecord, region: Tuple[int, int] = None, ax: mpl.axes = None
) -> plt.axes:
    """Plot Sanger chromatograph.

    region: include both start and end (1-based)
    """
    if seq is None:
        return ax

    if region is None:
        # turn into 0 based for better indexing
        region_start, region_end = 0, len(seq)
    else:
        region_start = max(region[0], 0)
        region_end = min(region[1], len(seq) - 1)

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(16, 6))

    _colors = defaultdict(
        lambda: "purple", {"A": "g", "C": "b", "G": "k", "T": "r"}
    )

    # Get signals
    peaks = seq.annotations["peak positions"]
    trace_x = seq.annotations["trace_x"]
    traces_y = [seq.annotations["channel " + str(i)] for i in range(1, 5)]
    bases = seq.annotations["channels"]

    xlim_left, xlim_right = peaks[region_start] - 1, peaks[region_end] + 0.5

    # subset peak and sequence
    # TODO: this might fix the bug
    peak_start = peaks[0]
    peak_zip = [
        (p, s)
        for i, (p, s) in enumerate(zip(peaks, seq))
        if region_start <= i <= region_end
    ]
    peaks, seq = list(zip(*peak_zip))

    # subset trace_x and traces_y together
    trace_zip = [
        (x + peak_start, *ys)
        for x, *ys in zip(trace_x, *traces_y)
        if xlim_left <= x <= xlim_right
    ]
    if not trace_zip:
        return ax
    trace_x, *traces_y = list(zip(*trace_zip))

    # Plot traces
    trmax = max(map(max, traces_y))
    for base in bases:
        trace_y = [1.0 * ci / trmax for ci in traces_y[bases.index(base)]]
        ax.plot(trace_x, trace_y, color=_colors[base], lw=2, label=base)
        ax.fill_between(
            trace_x, 0, trace_y, facecolor=_colors[base], alpha=0.125
        )

    # Plot bases at peak positions
    for i, peak in enumerate(peaks):
        #  LOGGER.debug(f"{i}, {peak}, {seq[i]}, {xlim_left + i}")
        ax.text(
            peak,
            -0.11,
            seq[i],
            color=_colors[seq[i]],
            va="center",
            ha="center",
            alpha=0.66,
            fontsize="x-large",
            fontweight="bold",
        )

    ax.set_ylim(bottom=-0.15, top=1.05)
    #  peaks[0] - max(2, 0.02 * (peaks[-1] - peaks[0])),
    #  right=peaks[-1] + max(2, 0.02 * (peaks[-1] - peaks[0])),
    ax.set_xlim(xlim_left + 0.5, xlim_right)
    ax.set_xticks(peaks)
    ax.set_xticklabels(list(range(region_start + 1, region_end + 2)))
    # hide y axis
    ax.set_yticklabels([])
    ax.get_yaxis().set_visible(False)
    # hide border
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # hide grid
    ax.grid(False)
    # set legend
    ax.legend(loc="upper left", bbox_to_anchor=(0.95, 0.99))
    return ax


def show_reference(
    query_record: SeqRecord,
    subject_record: SeqRecord,
    ax: mpl.axes,
    ref_central: Optional[int] = None,
) -> mpl.axes:
    """show the reference of the chromatograph.

    design: if location is not proviode, do the alignment first
    @param seq: input SeqRecord of ref
    """

    sitepairs = align_chromatograph(
        query_record, subject_record, ignore_ambig=True
    )
    sitepairs_indexing = {s.cf_pos: s for s in sitepairs}
    cf_sites = [int(i.get_text()) for i in ax.get_xticklabels()]
    matched_sitepairs = [sitepairs_indexing[pos] for pos in cf_sites]
    for i, peak in enumerate(ax.get_xticks()):
        ax.text(
            peak,
            1.05,
            matched_sitepairs[i].ref_base,
            color="dimgrey",
            va="bottom",
            ha="center",
            alpha=0.85,
            fontsize="xx-large",
            fontweight="bold",
            clip_on=False,
        )
        if ref_central is not None:
            ref_pos = matched_sitepairs[i].ref_pos - ref_central
        else:
            ref_pos = matched_sitepairs[i].ref_pos
        ax.text(
            peak,
            1.12,
            ref_pos,
            color="dimgrey",
            va="bottom",
            ha="center",
            alpha=0.85,
            fontsize="medium",
            fontweight="normal",
            clip_on=False,
        )
    return ax


def highlight_base(
    pos_highlight: int, seq: SeqRecord, ax: mpl.axes, passed_filter=True
) -> mpl.axes:
    """Highlight the area around a peak with a rectangle."""

    peaks = seq.annotations["peak positions"]
    peak = peaks[pos_highlight - 1]

    xmin, xmax = ax.get_xlim()
    if not xmin <= peak < xmax:
        raise ValueError("peak not within plot bounds")

    if pos_highlight == 1:
        xmin = -0.5
    else:
        xmin = 0.5 * (peaks[pos_highlight - 1] + peaks[pos_highlight - 2])

    if pos_highlight == len(peaks):
        xmax = -0.5
    else:
        xmax = 0.5 * (peaks[pos_highlight - 1] + peaks[pos_highlight])
    ymin, ymax = ax.get_ylim()

    if passed_filter:
        fcolor = "yellow"
    else:
        fcolor = "grey"
    rec = mpl.patches.Rectangle(
        (xmin, ymin),
        (xmax - xmin),
        (ymax - ymin),
        edgecolor="none",
        facecolor=fcolor,
        alpha=0.3,
    )
    ax.add_patch(rec)
    return ax


def annotate_mutation(mut: SitePair, seq: SeqRecord, ax) -> mpl.axes:
    """Annotate mutation pattern chromatograph position."""
    peaks = seq.annotations["peak positions"]
    peak = peaks[mut.cf_pos - 1]
    ax.text(
        peak,
        0.99,
        f"{mut.ref_base}{mut.ref_pos}{mut.cf_base}",
        color="c",
        fontsize="large",
        fontweight="bold",
        rotation=45,
        ha="center",
        va="center",
    )
    return ax
