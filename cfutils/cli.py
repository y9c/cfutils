#!/usr/bin/env python3
"""wrap cfutils into cli app"""

import click

from cfutils.run import do_mutation_calling


@click.group()
@click.option('--debug/--no-debug', default=False)
def cli(debug):
    """Chromatogram File Utils"""
    click.echo('Debug mode is %s' % ('on' if debug else 'off'))


@cli.command()
@click.option(
    '--query', prompt='QUERY (abi file): ', help='Query file in abi format')
@click.option(
    '--subject',
    prompt='SUBJECT (fasta file): ',
    help='Subject file in fasta format as ref')
@click.option(
    '--plot/--no-plot',
    default=False,
    help="Plot mutation location in chromatogram")
def mut(query, subject, plot):
    """do mutation calling"""
    do_mutation_calling(
        query_ab1_file=query,
        subject_fasta_file=subject,
        report_mut_info=True,
        report_mut_plot=plot,
        mut_info_file="./temp/test.tsv",
        mut_plot_file="./temp/test.pdf")


@cli.command()
def test():
    """test"""
    click.echo('Testing...')
