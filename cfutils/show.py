#!/usr/bin/env python3
'''
author:     Fabio Zanini
date:       09/12/13
content:    Plot functions for Sanger chromatographs.
modified:   By Ye Chang in 2018-05-14
'''

# Modules
import logging
import sys
from collections import defaultdict
from typing import List, Dict

from Bio.SeqRecord import SeqRecord

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError('cfutils requires Python 3.6+!')

LOGGER: logging.Logger = logging.getLogger()
HANDLER: logging.StreamHandler = logging.StreamHandler()
FORMATTER: logging.Formatter = logging.Formatter(
    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
HANDLER.setFormatter(FORMATTER)
LOGGER.addHandler(HANDLER)
LOGGER.setLevel(logging.DEBUG)
# LOGGER.setLevel(logging.INFO)

# Globals
BASES: List[str] = ['A', 'C', 'G', 'T']
COLORS = defaultdict(lambda: 'purple', {
    'A': 'g',
    'C': 'b',
    'G': 'k',
    'T': 'r'
})


# Functions
def plot_chromatograph(seq: SeqRecord, ax=None, xlim=None,
                       peaklim=None) -> None:
    '''Plot Sanger chromatograph'''

    if ax is None:
        import matplotlib.pyplot as plt
        plt.ion()
        fig, ax = plt.subplots(1, 1, figsize=(16, 6))

    if seq is None:
        ax.set_xlim(-2, 102)
        ax.set_ylim(-0.15, 1.05)
        return

    # Get signals
    traces = [seq.annotations['channel ' + str(i)] for i in range(1, 5)]
    peaks = seq.annotations['peak positions']
    BASES = seq.annotations['channels']
    x = seq.annotations['trace_x']

    # Limit to a region if necessary
    if (xlim is not None) or (peaklim is not None):
        if peaklim is not None:
            xlim = (peaks[min(len(peaks) - 1, peaklim[0])], peaks[min(
                len(peaks) - 1, peaklim[1] - 1)])
        ind = [(xi >= xlim[0]) and (xi <= xlim[1]) for xi in x]
        if not any(ind):
            return

        x = [xi for (indi, xi) in zip(ind, x) if indi]
        traces = [[ci for (indi, ci) in zip(ind, tr) if indi] for tr in traces]
        ind = [
            i for (i, ti) in enumerate(peaks)
            if (ti >= xlim[0]) and (ti <= xlim[1])
        ]
        peaks = peaks[ind[0]:ind[-1] + 1]
        seq = seq[ind[0]:ind[-1] + 1]

    # Plot traces
    trmax = max(map(max, traces))
    for base in BASES:
        y = [1.0 * ci / trmax for ci in traces[BASES.index(base)]]
        ax.plot(x, y, color=COLORS[base], lw=2, label=base)
        ax.fill_between(x, 0, y, facecolor=COLORS[base], alpha=0.125)

    # Plot bases at peak positions
    LOGGER.debug(seq)
    for i, peak in enumerate(peaks):
        LOGGER.debug("%i, %f, %s, %i", i, peak, seq[i], xlim[0] + i)
        ax.text(
            peak,
            -0.11,
            seq[i],
            color=COLORS[seq[i]],
            horizontalalignment='center')

    ax.set_ylim(ymin=-0.15, ymax=1.05)
    ax.set_xlim(
        xmin=peaks[0] - max(2, 0.02 * (peaks[-1] - peaks[0])),
        xmax=peaks[-1] + max(2, 0.02 * (peaks[-1] - peaks[0])))
    ax.set_yticklabels([])
    ax.set_xticks(peaks)
    ax.set_xticklabels(list(range(xlim[0] - 2, xlim[0] + len(peaks) - 2)))
    ax.grid(False)
    ax.legend(loc='upper left', bbox_to_anchor=(0.93, 0.99))


def closest_peak(pos_highlight: int, seq: SeqRecord) -> Dict:
    peaks = seq.annotations['peak positions']
    (i, peak) = min(
        enumerate(peaks), key=lambda x: abs(x[0] + 1 - pos_highlight))
    return {'index': i, 'peak': peak}


def peak_position(i, seq):
    return seq.annotations['peak positions'][i]


def highlight_base(pos_highlight, seq, ax):
    '''Highlight the area around a peak with a rectangle'''

    #  trace = seq.annotations['channel 1']
    peaks = seq.annotations['peak positions']

    peak_obj = closest_peak(pos_highlight, seq)
    i = peak_obj['index']
    peak = peak_obj['peak']

    xmin, xmax = ax.get_xlim()
    if not xmin <= peak < xmax:
        raise ValueError('peak not within plot bounds')

    if i == 0:
        xmin = -0.5
    else:
        xmin = 0.5 * (peaks[i - 1] + peak)

    if i == len(peaks) - 1:
        xmax = peak + 0.5
    else:
        xmax = 0.5 * (peak + peaks[i + 1])

    ymin, ymax = ax.get_ylim()

    from matplotlib.patches import Rectangle
    rec = Rectangle(
        (xmin, ymin), (xmax - xmin), (ymax - ymin),
        edgecolor='none',
        facecolor='yellow',
        alpha=0.3)
    ax.add_patch(rec)
    return {'index': i, 'peak': peak, 'rec': rec}
