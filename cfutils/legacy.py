#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
#
# Distributed under terms of the MIT license.

"""Chromatogram File Utils.

- update in 20190405
"""

import os
import re
import subprocess
import sys
from io import StringIO
from optparse import OptionParser
from subprocess import PIPE, STDOUT, Popen

from Bio import SeqIO, pairwise2
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA


# get option
def check_options(options):
    if (
        options.subject_fasta
        and options.file_name_raw
        and options.company_format
    ):
        return True
    return False


def get_options():
    o = OptionParser()
    o.add_option(
        "--raw",
        dest="file_name_raw",
        help="*Complessing Raw sequencing result file",
    )
    o.add_option(
        "--map",
        dest="corresponding_table",
        help="table file tha Map sequence, primer and reference",
    )
    o.add_option(
        "--company",
        dest="company_format",
        type="choice",
        default="invitrogen",
        choices=["invitrogen", "bgi", "generay", "genewiz"],
    )
    o.add_option(
        "--subject",
        dest="subject_fasta",
        help="*Subject fasta file that has only 1 sequence to align all the query sequences against. Most likely a reference gene fasta file.",
    )
    o.add_option(
        "--skip-ambiguous",
        dest="ignore_ambig",
        action="store_false",
        default=True,
        help="Specify this option to skip ambiguous bases so they are not counted",
    )
    options, args = o.parse_args()

    if not check_options(options):
        o.print_help()
        sys.exit(-1)

    return options


# convert the complessing file into fasta file and abi folder


def decomplessing_file(file_name_raw):
    print(file_name_raw)
    folder_name = ".".join(file_name_raw.split(".")[0:-1])
    print(folder_name)
    cp1 = subprocess.Popen("rm -rf '{}'".format(folder_name), shell=True)
    cp1.wait()
    print("rm -rf '{}'".format(folder_name))
    cp2 = subprocess.Popen("mkdir '{}'".format(folder_name), shell=True)
    cp2.wait()
    if file_name_raw[-3:] == "rar":
        cp3 = subprocess.Popen(
            "rar x '{}' '{}'".format(file_name_raw, folder_name), shell=True
        )
    if file_name_raw[-3:] == "zip":
        cp3 = subprocess.Popen(
            "unzip -q -n '{}' -d '{}'".format(file_name_raw, folder_name),
            shell=True,
        )
        print("unzip -q -n '{}' -d '{}'".format(file_name_raw, folder_name))
    if file_name_raw[-6:] == "tar.gz":
        cp3 = subprocess.Popen(
            "tar -zcxf '{}' -C '{}'".format(file_name_raw, folder_name),
            shell=True,
        )
    cp3.wait()
    return folder_name


######################
# parse files method #
######################


def parse_ab1(ab1_file):
    with open(ab1_file, "rb") as handle:
        seq = next(SeqIO.parse(handle, "abi"))
    return seq


def parse_seq(seq_file):
    handle = open(seq_file, "r")
    return handle.read().replace("^M$", "")
    handle.close()


##################
# classify files #
##################


def classify(folder_name, company_format, corresponding_table):
    map_table = {}
    for l in open(corresponding_table, "r").readlines()[1:]:
        ll = re.split("\W+", l)
        map_table[ll[0]] = ll[1]

    group_dic_of_file_path = {}
    for (path, subdirs, files) in os.walk(folder_name):
        for name in files:
            file_path = os.path.join(path, name)
            if (
                name[-4:] == ".abl"
                or name[-4:] == ".ab1"
                or name[-4:] == ".abi"
            ):
                if company_format == "invitrogen" and len(name.split(".")) > 5:
                    ORDERID, RUNID, SAMPLENAME, PRIMERNAME, REACTIONID, FILETYPE = name.split(
                        "."
                    )
                    ref_name = map_table[SAMPLENAME]
                    if ref_name in group_dic_of_file_path:
                        group_dic_of_file_path[ref_name].append(file_path)
                    else:
                        group_dic_of_file_path[ref_name] = [file_path]
    return group_dic_of_file_path


#########################
# pairwise align method #
#########################


def pairwise_align(seq1, seq2, ignore_ambig):
    alignments = pairwise2.align.globalms(seq1, seq2, 1, -2, -5, -2)
    get_muts(alignments[0][0], alignments[0][1], ignore_ambig)


def tcoffee_align(seq1, seq2, ignore_ambig):
    seqs = StringIO()
    seqs.write(seq1.format("fasta"))
    seqs.write(seq2.format("fasta"))
    cmd = "t_coffee -in stdin -output fasta -outfile stdout -quiet"
    p = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    stdeo, stdin = p.communicate(input=seqs.getvalue())
    s = StringIO(stdeo)

    try:
        p = SeqIO.parse(s, "fasta")
        ref = next(p)
        s = next(p)
    except Exception:
        return (stdeo, -1)

    mutations = get_muts(ref.seq, s.seq, ignore_ambig)

    return (mutations, 1)


def is_ambig(base):
    """If a given base is ambiguous or not."""
    return base.upper() not in IUPACUnambiguousDNA.letters


def get_muts(sub_align, query_align, ignore_ambig):
    mutations = []
    count = 1
    for q, s in zip(query_align, sub_align):
        tmp = "Q: %s S: %s Pos: %s" % (q, s, count)

        if not ignore_ambig and (is_ambig(q) or is_ambig(s)):
            sys.stderr.write("Skipping ambiguous base\n")
            sys.stderr.write(tmp + "\n")
        elif q != s:
            mutations.append(tmp)
        count += 1
    return mutations


def fasta_into_dic(fasta_file):
    fasta_dic = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dic[record.name] = record
    return fasta_dic


##################
# deal each file #
##################


def align(
    file_name_raw,
    subject_fasta,
    company_format,
    corresponding_table,
    ignore_ambig,
):
    folder_name = decomplessing_file(file_name_raw)
    group_dic_of_file_path = classify(
        folder_name, company_format, corresponding_table
    )
    ref_fasta_dic = fasta_into_dic(subject_fasta)
    count = 1
    for ref_name in group_dic_of_file_path:
        ref_sequence = ref_fasta_dic[ref_name]
        for file_path in group_dic_of_file_path[ref_name]:
            query_seq_abi = parse_ab1(file_path)
            query_seq = query_seq_abi

            print(
                "the {}th file :----------->\n{}& {}".format(
                    count, ref_name, file_path
                )
            )
            count += 1
            mutations = []
            mutations, r1 = tcoffee_align(
                ref_sequence, query_seq, ignore_ambig
            )
            for m in mutations:
                print(m)


# main part of the program
def main(options):
    align(
        options.file_name_raw,
        options.subject_fasta,
        options.company_format,
        options.corresponding_table,
        options.ignore_ambig,
    )


if __name__ == "__main__":
    options = get_options()
    main(options)
