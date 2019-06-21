#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-05-27 20:19


"""align two sequence with ref by blast.

Use 1-based for all the position
"""

import shutil
import sys
import tempfile
from dataclasses import dataclass
from io import StringIO
from subprocess import PIPE, STDOUT, Popen
from typing import List, Optional, Tuple

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA, ambiguous_dna
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blast
from Bio.Blast.Record import HSP, Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .utils import get_logger

LOGGER = get_logger(__name__)


@dataclass
class SitePair:
    """Object for storing align pair at mutation site."""

    ref_pos: int
    ref_base: str
    cf_pos: int
    cf_base: str
    qual_site: Optional[int] = None
    qual_local: Optional[int] = None

    def __repr__(self):
        return (
            f"{self.ref_base}({self.ref_pos})->{self.cf_base}({self.cf_pos})"
        )


def is_ambig(base):
    """If a given base is ambiguous or not."""
    return base.upper() not in IUPACUnambiguousDNA.letters + "-"


def rc_seq(seq):
    """reverse_complement sequence in str."""
    return str(Seq(seq, ambiguous_dna).reverse_complement())


def which_blast():
    """check blastn exist."""
    return shutil.which("blastn")


def run_blast(
    query_record: SeqRecord,
    subject_record: SeqRecord,
    ignore_ambig: bool = False,
) -> List[SitePair]:
    """run blastn."""
    if not which_blast():
        LOGGER.error("Please ensure blastn is installed and in your PATH")
        sys.exit(1)
    with tempfile.NamedTemporaryFile() as query_file, tempfile.NamedTemporaryFile() as subject_file:
        query_file.write(query_record.format("fasta").encode())
        query_file.seek(0)
        subject_file.write(subject_record.format("fasta").encode())
        subject_file.seek(0)

        cline = blast(
            cmd="blastn",
            query=query_file.name,
            subject=subject_file.name,
            #  gapopen=5,
            #  gapextend=2,
            gapopen=10,
            gapextend=4,
            reward=3,
            penalty=-4,
            outfmt=5,
        )
        LOGGER.debug(f"Run blastn with this command:\n{cline}")
        with Popen(str(cline).split(), stdout=PIPE, stderr=STDOUT) as proc:
            stdeo, _ = proc.communicate()
            #  LOGGER.debug(stdeo.decode())

    # NOTE: this SitePair object is without qual
    hsp = parse_blast(stdeo, ignore_ambig=ignore_ambig)
    return hsp


def parse_blast(output, ignore_ambig=False):
    """parse blastn output.

    @return hsp: HSP (high-scoring pair)
    """
    blast_output = StringIO(output.decode())

    try:
        blast_records: Alignment = NCBIXML.read(blast_output)
    except ValueError as err:
        if (
            blast_output.getvalue()
            == "BLAST engine error: XML formatting is only supported for a database search"
        ):
            LOGGER.warning(
                "Please ensure that you are using the latest version of blastn"
                "You may need to update your environment's PATH variable"
            )
        raise err

    # TODO: check and score alignment result
    # the second hit?
    # the alignment score
    if not blast_records.alignments:
        LOGGER.info(
            f"Can not find alignment of the ab1 file with the ref sequence!"
        )
        return []
    alignment = blast_records.alignments[0]
    # only return the first result
    hsp = alignment.hsps[0]
    return hsp


def get_pairs(
    hsp: HSP, ignore_ambig: bool = False, mut_only: bool = False
) -> List[SitePair]:
    """get_pairs for parse SitePair object from hsp object.

    # can be used as `get_muts`, when mut_only == True

    NOTE: HSP object is 1-based
    @param: hsp:
    @param: ignore_ambig:
    @rtype: List[List[Union[int, str]]]
    """
    sitepairs = []
    if hsp.sbjct_start < hsp.sbjct_end:
        sbjct_loc = hsp.sbjct_start
        query_loc = hsp.query_start
        for query_base, sbjct_base in zip(hsp.query, hsp.sbjct):

            if ignore_ambig and is_ambig(sbjct_base):
                LOGGER.debug(f"Skipping ambiguous base: {sbjct_base}")
            elif mut_only and query_base == sbjct_base:
                LOGGER.debug(f"Skipping non mut base: {sbjct_base}")
            else:
                sitepairs.append(
                    SitePair(sbjct_loc, sbjct_base, query_loc, query_base)
                )
            if sbjct_base != "-":
                sbjct_loc += 1
            if query_base != "-":
                query_loc += 1
    else:
        sbjct_loc = hsp.sbjct_end
        query_loc = hsp.query_end
        for query_base, sbjct_base in zip(
            rc_seq(hsp.query), rc_seq(hsp.sbjct)
        ):

            if ignore_ambig and is_ambig(sbjct_base):
                LOGGER.info("Skipping ambiguous base\n")
            elif mut_only and query_base == sbjct_base:
                LOGGER.debug(f"Skipping non mut base: {sbjct_base}")
            else:
                sitepairs.append(
                    SitePair(sbjct_loc, sbjct_base, query_loc, query_base)
                )
            if sbjct_base != "-":
                sbjct_loc += 1
            if query_base != "-":
                query_loc -= 1
    return sitepairs


def get_quality(
    pos: int, query_record: SeqRecord, flank_base_num=0
) -> Tuple[int, int]:
    """get quality of site and local region.

    change flank_base_num to number gt 0 to get mean qual within region
    """
    qual = query_record.letter_annotations["phred_quality"]
    qual_site = qual[pos - 1]
    qual_flank = qual[
        max(0, pos - 1 - flank_base_num) : min(len(qual), pos + flank_base_num)
    ]
    qual_local = int(sum(qual_flank) / len(qual_flank))
    return qual_site, qual_local


def align(
    query_record: SeqRecord, subject_record: SeqRecord, ignore_ambig=False
) -> List[SitePair]:
    """run align.

    @return: list of SitePair about all sites
    """
    hsp = run_blast(query_record, subject_record, ignore_ambig=ignore_ambig)
    sitepairs = get_pairs(hsp, ignore_ambig=ignore_ambig, mut_only=False)
    LOGGER.info(
        f"{query_record.description}: Total aligned number: {len(sitepairs)}"
    )
    for site in sitepairs:
        site.qual_site, site.qual_local = get_quality(
            site.cf_pos, query_record, flank_base_num=5
        )
        LOGGER.debug(f"{site}\tlocal:{site.qual_local}\tsite:{site.qual_site}")
    return sitepairs


def call_mutations(
    query_record: SeqRecord, subject_record: SeqRecord, ignore_ambig=False
) -> List[SitePair]:
    """run align and call mutations.

    @return: list of SitePair about mutation sites
    """
    hsp = run_blast(query_record, subject_record, ignore_ambig=ignore_ambig)
    mutations = get_pairs(hsp, ignore_ambig=ignore_ambig, mut_only=True)
    LOGGER.info(
        f"{query_record.description}: Total mutation number: {len(mutations)}"
    )
    for site in mutations:
        site.qual_site, site.qual_local = get_quality(
            site.cf_pos, query_record, flank_base_num=5
        )
        LOGGER.debug(f"{site}\tlocal:{site.qual_local}\tsite:{site.qual_site}")
    return mutations
