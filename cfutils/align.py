#!/usr/bin/env python3
"""
align two sequence with ref by blast
"""

import logging
import sys
from io import StringIO
from subprocess import PIPE, STDOUT, Popen

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA, ambiguous_dna
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blast
from Bio.Seq import Seq

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
# LOGGER.setLevel(logging.DEBUG)
LOGGER.setLevel(logging.INFO)


def run_blast(seq, subject_fasta, ignore_ambig=True):
    cline = blast(
        cmd='./bin/blastn',
        subject=subject_fasta,
        gapopen=5,
        gapextend=2,
        reward=3,
        penalty=-4,
        outfmt=5)
    #, out='blastoutput.xml' ) #outfmt="6 qseqid sseqid evalue slen mismatch" )#, out='wrair2368t_pb2.xml' )
    try:
        p = Popen(str(cline).split(), stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    except OSError:
        LOGGER.error("Please ensure blastn is installed and in your PATH")
        sys.exit(1)
    stdeo, stdin = p.communicate(input=seq.format('fasta').encode())
    LOGGER.debug(stdeo.decode())

    return parse_blast(seq, stdeo, ignore_ambig=ignore_ambig)


def parse_blast(seq, output, ignore_ambig=True):
    blast_output = StringIO(output.decode())

    try:
        blast_records = NCBIXML.read(blast_output)
    except ValueError as e:
        LOGGER.info("-----Blast output------")
        LOGGER.info(blast_output.getvalue())
        if blast_output.getvalue(
        ) == "BLAST engine error: XML formatting is only supported for a database search":
            LOGGER.warning(
                "Please ensure that you are using the latest blastx version of blastn"
            )
            LOGGER.warning(
                "You may need to update your environment's PATH variable")
        raise e

    try:
        alignment = blast_records.alignments[0]
    except:
        return (output, -1)

    hsp = alignment.hsps[0]

    #  mutations = get_muts(hsp.query, hsp.sbjct, ignore_ambig=ignore_ambig)
    mutations = get_muts(hsp, ignore_ambig=ignore_ambig)

    return (mutations, 1)


def is_ambig(base):
    """
        If a given base is ambiguous or not
    """
    return base.upper() not in IUPACUnambiguousDNA.letters


def rc_seq(seq):
    """
        reverse_complement sequence in str
    """
    return str(Seq(seq, ambiguous_dna).reverse_complement())


def get_muts(hsp, ignore_ambig=True):
    mutations = []
    if hsp.sbjct_start < hsp.sbjct_end:
        sbjct_loc = hsp.sbjct_start
        query_loc = hsp.query_start
        for q, s in zip(hsp.query, hsp.sbjct):
            mutation_symbol = f"RefLocation: {sbjct_loc}\tRefBase: {s}\tCfLocation: {query_loc}\tMutBase: {q}"

            if not ignore_ambig and (is_ambig(q) or is_ambig(s)):
                LOGGER.info("Skipping ambiguous base\n")
                LOGGER.info(mutation_symbol)
            elif q != s:
                mutations.append(mutation_symbol)
            if s != "-":
                sbjct_loc += 1
            if q != "-":
                query_loc += 1
    else:
        sbjct_loc = hsp.sbjct_end
        query_loc = hsp.query_end
        for q, s in zip(rc_seq(hsp.query), rc_seq(hsp.sbjct)):
            mutation_symbol = f"RefLocation: {sbjct_loc}\tRefBase: {s}\tCfLocation: {query_loc}\tMutBase: {q}"

            if not ignore_ambig and (is_ambig(q) or is_ambig(s)):
                LOGGER.info("Skipping ambiguous base\n")
                LOGGER.info(mutation_symbol)
            elif q != s:
                mutations.append(mutation_symbol)
            if s != "-":
                sbjct_loc += 1
            if q != "-":
                query_loc -= 1
    return mutations


def align(query_record, subject_fasta, ignore_ambig=True):
    if not ignore_ambig:
        LOGGER.info("Ignoring ambiguous bases")

    mutations = []
    mutations, status_blastn = run_blast(
        query_record, subject_fasta, ignore_ambig=ignore_ambig)
    LOGGER.info("%s: Total mutations: %s" % (query_record.description,
                                             len(mutations)))
    for m in mutations:
        print(m)
