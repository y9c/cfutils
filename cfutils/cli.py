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
def mut(query, subject):
    """do mutation calling"""
    do_mutation_calling(query, subject)


@cli.command()
def test():
    """test"""
    click.echo('Testing...')
