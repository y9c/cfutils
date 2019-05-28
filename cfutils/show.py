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
from typing import List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord

from .align import SitePair
from .utils import get_logger

LOGGER = get_logger(__name__)


_BASES = ["A", "C", "G", "T"]
_COLORS = defaultdict(
    lambda: "purple", {"A": "g", "C": "b", "G": "k", "T": "r"}
)


def plot_chromatograph(seq: SeqRecord, ax=None, region: Tuple = None) -> None:
    """Plot Sanger chromatograph.

    region: include both start and end
    """

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(16, 6))

    if seq is None:
        ax.set_xlim(-2, 102)
        ax.set_ylim(-0.15, 1.05)
        return

    # Get signals
    traces = [seq.annotations["channel " + str(i)] for i in range(1, 5)]
    peaks = seq.annotations["peak positions"]
    BASES = seq.annotations["channels"]
    x = seq.annotations["trace_x"]

    # Limit to a region if necessary
    if region is not None:
        xlim = (
            peaks[min(len(peaks), region[0]) - 1],
            peaks[min(len(peaks), region[1]) - 1],
        )
        ind = [(xi >= xlim[0]) and (xi <= xlim[1]) for xi in x]
        if not any(ind):
            return

        x = [xi for (indi, xi) in zip(ind, x) if indi]
        traces = [[ci for (indi, ci) in zip(ind, tr) if indi] for tr in traces]
        ind = [
            i
            for (i, ti) in enumerate(peaks)
            if (ti >= xlim[0]) and (ti <= xlim[1])
        ]
        peaks = peaks[ind[0] : ind[-1] + 1]
        seq = seq[ind[0] : ind[-1] + 1]

    # Plot traces
    trmax = max(map(max, traces))
    for base in BASES:
        y = [1.0 * ci / trmax for ci in traces[BASES.index(base)]]
        ax.plot(x, y, color=_COLORS[base], lw=2, label=base)
        ax.fill_between(x, 0, y, facecolor=_COLORS[base], alpha=0.125)

    # Plot bases at peak positions
    LOGGER.debug(seq)
    for i, peak in enumerate(peaks):
        LOGGER.debug("%i, %f, %s, %i", i, peak, seq[i], xlim[0] + i)
        ax.text(
            peak,
            -0.11,
            seq[i],
            color=_COLORS[seq[i]],
            horizontalalignment="center",
        )

    ax.set_ylim(bottom=-0.15, top=1.05)
    ax.set_xlim(
        left=peaks[0] - max(2, 0.02 * (peaks[-1] - peaks[0])),
        right=peaks[-1] + max(2, 0.02 * (peaks[-1] - peaks[0])),
    )
    ax.set_yticklabels([])
    ax.set_xticks(peaks)
    ax.set_xticklabels(list(range(region[0], region[0] + len(peaks))))
    # hide border
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # hide y axis
    ax.get_yaxis().set_visible(False)
    ax.grid(False)
    ax.legend(loc="upper left", bbox_to_anchor=(0.95, 0.99))


def highlight_base(
    pos_highlight: int, seq: SeqRecord, ax, passed_filter=True
) -> Tuple:
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
    return (xmin, xmax)


def annotate_mutation(mut: SitePair, seq: SeqRecord, ax) -> None:
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
