#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Chromatogram File Utils.

wrap cfutils into cli app
- update in 20190405
"""

import click

from cfutils.run import report_mutation


@click.group()
@click.option("--debug/--no-debug", default=False)
def cli(debug):
    """Chromatogram File Utils."""
    if debug:
        click.echo("Debug mode is on")


# call mutation
@cli.command()
@click.option(
    "--query", prompt="QUERY (abi file): ", help="Query file in abi format"
)
@click.option(
    "--subject",
    prompt="SUBJECT (fasta file): ",
    help="Subject file in fasta format as ref",
)
@click.option(
    "--outdir", default=None, required=False, help="Output directory"
)
@click.option(
    "--outbase", default=None, required=False, help="Output basename"
)
@click.option(
    "--aligned/--mutated",
    default=False,
    help="Report all aligned sites or mutation sites only",
)
@click.option(
    "--plot/--no-plot",
    default=False,
    help="Generate figure of mutation in chromatogram.",
)
def mut(query, subject, outdir, outbase, aligned, plot):
    """do mutation calling, then report in tsv and pdf."""
    report_mutation(
        query_ab1_file=query,
        subject_fasta_file=subject,
        output_dir=outdir,
        file_basename=outbase,
        report_all_sites=aligned,
        report_mut_plot=plot,
    )


# test
@cli.command()
def test():
    """test."""
    click.echo("Testing...")
