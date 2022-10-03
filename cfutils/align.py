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

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blast
from Bio.Blast.Record import HSP, Alignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from parasail import blosum62, sw_trace_scan_sat

from .utils import chunked_lines, get_logger

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


try:
    from parasail import (
        blosum62,
        nw_banded,
        nw_trace_scan_sat,
        ssw,
        sw_trace_scan_sat,
    )
except ImportError:
    import warnings

    warnings.warn(
        "ERROR!!! could not import parasail... will not be able to align"
    )


def align_seqs(query, target, gop=5, gep=10):
    """basic parasail global alignment of two sequences result is wrapped in
    ParasailAlignment Class."""
    query = str(query)
    target = str(target)
    result = sw_trace_scan_sat(target, query, gop, gep, blosum62)
    return ParasailAlignment(result)


def parental_numbering(aseq1, aseq2):
    """given two ALIGNED sequences, return a 'position list' for the second
    sequence based on the parental sequence."""
    idx = 1
    numlist = []
    insertchars = "abcdefghijklmnopqrstuvwxyz"
    insertidx = 0
    for s1, s2 in zip(aseq1, aseq2):
        if s2 == "-":
            idx += 1
            continue
        if s1 == "-":
            numlist.append(
                str(idx - 1) + insertchars[insertidx % len(insertchars)]
            )
            insertidx += 1
            continue
        insertidx = 0
        numlist.append(str(idx))
        idx += 1
    return numlist


class ParasailAlignment:
    """Convenience class to wrap the results of a parasail alignment."""

    def __init__(self, result):
        self.cigar = result.cigar.decode
        if isinstance(self.cigar, bytes):
            self.cigar = self.cigar.decode()
        # confusing nomenclature is for consistency with scikit-bio
        # where "query" is the initial sequence
        self.target = result.query
        self.query = result.ref
        self.score = result.score
        self.start_target = 0
        self.start_query = 0
        self.end_target = result.end_query
        self.end_query = result.end_ref
        self._mutations = None

    def _tuples_from_cigar(self):
        tuples = []
        length_stack = []
        for character in self.cigar:
            if character.isdigit():
                length_stack.append(character)
            else:
                tuples.append((int("".join(length_stack)), character))
                length_stack = []
        return tuples

    @property
    def cigar_tuple(self):
        if hasattr(self, "_cigar_tuple"):
            return self._cigar_tuple
        self._cigar_tuple = self._tuples_from_cigar()
        return self._cigar_tuple

    def __repr__(self):
        return str(self)

    def __str__(self):
        a = chunked_lines(self.aligned_target_sequence(), spacer="")
        b = chunked_lines(self.aligned_query_sequence(), spacer="")
        out = []
        for t, q in zip(a, b):
            out.append(q)
            out.append(
                "".join(
                    [
                        "*" if x != y else (" " if x == " " else "|")
                        for x, y in zip(t, q)
                    ]
                )
            )
            out.append(t + "\n")
        return "\n".join(out)

    def __iter__(self):
        yield self.aligned_query_sequence()
        yield self.aligned_target_sequence()

    def as_mutations(self, reference=None):
        from .mutations import MutationSet, _get_aligned_muts

        seq1, seq2 = self
        mutstring = "/".join(_get_aligned_muts(seq1, seq2))
        if reference is not None:
            return MutationSet(
                mutstring, parental_numbering(*align_seqs(reference, seq1))
            )
        return MutationSet(mutstring)

    def print_alignment(self, max_length=80):
        print(
            self.aligned_query_sequence()
            + "\n"
            + self.aligned_target_sequence()
        )

    def aligned_query_sequence(self):
        return self._get_aligned_sequence(self.query, "I")

    def aligned_target_sequence(self):
        return self._get_aligned_sequence(self.target, "D")

    def _get_aligned_sequence(self, seq, gap_type, gap_char="-", eq_char="="):
        # assume zero based
        # gap_type is 'D' when returning aligned query sequence
        # gap_type is 'I' when returning aligned target sequence
        aligned_sequence = ""
        index = 0
        for length, symbol in self.cigar_tuple:
            if symbol in (eq_char, "X"):
                end = length + index
                aligned_sequence += seq[index:end]
                index += length
            elif symbol == gap_type:
                aligned_sequence += gap_char * length
            elif symbol in ("D", "I"):
                end = length + index
                aligned_sequence += seq[index:end]
                index += length
        return aligned_sequence

    def _get_aligned_tuple(self):
        print(self.query)
        print(self.target)
        print(self.cigar_tuple)
        # assume zero based
        # gap_type is 'D' when returning aligned query sequence
        # gap_type is 'I' when returning aligned target sequence
        aligned_sequence = ""
        index = 0
        for length, symbol in self.cigar_tuple:
            if symbol in ("=", "X"):
                end = length + index
                aligned_sequence += seq[index:end]
                index += length
            elif symbol == "D":
                aligned_sequence += gap_char * length
            elif symbol == "I":
                end = length + index
                aligned_sequence += seq[index:end]
                index += length
        return aligned_sequence

    @classmethod
    def from_seqs(cls, query, target, **kwargs):
        return align_seqs(query, target, **kwargs)


def is_ambig(base):
    """If a given base is ambiguous or not."""
    return base.upper() not in "ATGC-"


def rc_seq(seq):
    """reverse_complement sequence in str."""
    return str(Seq(seq).reverse_complement())


def which_blast():
    """check blastn exist."""
    return shutil.which("blastn")


def parse_blast(output):
    """parse blastn output.

    @return hsp: HSP (high-scoring pair)
    """
    blast_output = StringIO(output.decode())

    try:
        blast_records: Alignment = NCBIXML.read(blast_output)
    except ValueError as err:
        if (
            blast_output.getvalue() == "BLAST engine error: "
            "XML formatting is only supported for a database search"
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
        #  exit if not aligned
        #  return []
        sys.exit(0)
    alignment = blast_records.alignments[0]
    # only return the first result
    hsp = alignment.hsps[0]
    if hsp.align_length < 30:
        LOGGER.info(
            f"alignment of the ab1 file with the ref sequence is < 30bp!"
        )
        #  exit if alignment is not short
        #  return []
        sys.exit(1)
    LOGGER.info(
        "The length and span of alignment: "
        f"{hsp.align_length} ({hsp.sbjct_start}-{hsp.sbjct_end})"
    )
    return hsp


def run_blast(query_record: SeqRecord, subject_record: SeqRecord) -> HSP:
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

    # NOTE: this SitePair object is without qual in this step
    hsp = parse_blast(stdeo)
    return hsp


def get_pairs(hsp: HSP, ignore_ambig: bool = False) -> List[SitePair]:
    """get_pairs for parse SitePair object from hsp object.

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


def align_chromatograph(
    query_record: SeqRecord, subject_record: SeqRecord, ignore_ambig=False
) -> List[SitePair]:
    """run align.

    @return: list of SitePair about all sites
    """
    hsp = run_blast(query_record, subject_record)
    sitepairs = get_pairs(hsp, ignore_ambig=ignore_ambig)
    LOGGER.info(f"{query_record.name}: Total aligned number: {len(sitepairs)}")
    for site in sitepairs:
        site.qual_site, site.qual_local = get_quality(
            site.cf_pos, query_record, flank_base_num=5
        )
        LOGGER.debug(f"{site}\tlocal:{site.qual_local}\tsite:{site.qual_site}")
    return sitepairs


def call_mutations(
    query_record: SeqRecord,
    subject_record: SeqRecord,
    ignore_ambig: bool = False,
    report_all_sites: bool = False,
) -> List[SitePair]:
    """run align and call mutations.

    @return: list of SitePair about mutation sites
    """
    print(".....\n\n")
    print(query_record.seq, subject_record)
    aa = align_seqs(query_record.seq, subject_record.seq)
    print(aa)
    aa._get_aligned_tuple()
    sitepairs = align_chromatograph(
        query_record, subject_record, ignore_ambig=ignore_ambig
    )
    mutations = []
    for site in sitepairs:
        if report_all_sites:
            mutations.append(site)
            LOGGER.debug(f"Site ({site}) is reported!")
        else:
            if site.ref_base != site.cf_base:
                mutations.append(site)
                LOGGER.debug(f"Site ({site}) is with mutation!")
    if not report_all_sites:
        LOGGER.info(
            f"{query_record.name}: Total mutation number: {len(mutations)}"
        )
    return mutations
