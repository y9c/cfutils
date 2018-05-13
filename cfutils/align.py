# Biopython imports

from Bio.Blast.Applications import NcbiblastnCommandline as blast
from Bio.Blast import NCBIXML
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA

from StringIO import StringIO
from subprocess import Popen, PIPE, STDOUT
import sys

# Wether or not ambiguous bases should be called different or not
ignore_ambig = True

def run_blast( seq, subject_fasta ):
    cline = blast( cmd='blastn', subject=subject_fasta, gapopen=5, gapextend=2, reward=1, penalty=-2, outfmt=5 )#, out='blastoutput.xml' ) #outfmt="6 qseqid sseqid evalue slen mismatch" )#, out='wrair2368t_pb2.xml' )
    try:
        p = Popen( str(cline).split(), stdout=PIPE, stdin=PIPE, stderr=STDOUT )
    except OSError:
        print "Please ensure blastn is installed and in your PATH"
        sys.exit( 1 )
    stdeo, stdin = p.communicate( input=seq.format( 'fasta' ) )

    return parse_blast( seq, stdeo )

def parse_blast( seq, output ):
    blast_output = StringIO( output )

    try:
        blast_records = NCBIXML.read( blast_output )
    except ValueError as e:
        sys.stderr.write( "-----Blast output------" )
        sys.stderr.write( blast_output.getvalue() )
        if blast_output.getvalue() == "BLAST engine error: XML formatting is only supported for a database search":
            sys.stderr.write( "Please ensure that you are using the latest blastx version of blastn" )
            sys.stderr.write( "You may need to update your environment's PATH variable" )
        raise e

    try:
        alignment = blast_records.alignments[0]
    except:
        return (output, -1)

    hsp = alignment.hsps[0]

    mutations = get_muts( hsp.query, hsp.sbjct )

    return (mutations,1)

def pairwise_align( seq1, seq2 ):
    alignments = pairwise2.align.globalms( seq1, seq2, 1, -2, -5, -2 )
    count = 1
    get_muts( alignments[0][0], alignments[0][1] )

def tcoffee_align( seq1, seq2 ):
    seqs = StringIO()
    seqs.write( seq1.format( 'fasta' ) )
    seqs.write( seq2.format( 'fasta' ) )
    cmd = "t_coffee -in stdin -output fasta -outfile stdout -quiet"
    p = Popen( cmd.split(), stdout=PIPE, stdin=PIPE, stderr=STDOUT )
    stdeo, stdin = p.communicate( input = seqs.getvalue() )
    s = StringIO( stdeo )

    try:
        p = SeqIO.parse( s, 'fasta' )
        ref = p.next()
        s = p.next()
    except Exception:
        return (stdeo, -1)

    mutations = get_muts( ref.seq, s.seq )

    return (mutations,1)

def is_ambig( base ):
    """
        If a given base is ambiguous or not
    """
    return base.upper() not in IUPACUnambiguousDNA.letters

def get_muts( sub_align, query_align ):
    mutations = []
    count = 1
    for q,s in zip( query_align, sub_align ):
        tmp = "Q: %s S: %s Pos: %s" % (q, s, count)

        if not ignore_ambig and (is_ambig( q ) or is_ambig( s )):
            sys.stderr.write( "Skipping ambiguous base\n" )
            sys.stderr.write( tmp + "\n" )
        elif q != s:
            mutations.append( tmp )
        count += 1
    return mutations

def count_seqs( file_name ):
    total_seq = 0
    fh = open( file_name )
    for l in fh:
        if l[0] == '>':
            total_seq += 1
    fh.close()

    return total_seq

def align( query_fasta, subject_fasta, ignore_ambiguous ):
    global ignore_ambig
    ignore_ambig = ignore_ambiguous
    if not ignore_ambig:
        sys.stderr.write( "Ignoring ambiguous bases" )

    refp = SeqIO.parse( subject_fasta, 'fasta' )
    ref = refp.next()

    fhq = SeqIO.parse( query_fasta, 'fasta' )

    total_seq = count_seqs( query_fasta )

    tcount = 1
    for seq in fhq:
        mutations = []
        sys.stderr.write( "Gathering mutations for sequence %s of %s\n" % (tcount, total_seq) )

        mutations, r1 = run_blast( seq, subject_fasta )
        if r1 == -1:
            sys.stderr.write( "%s failed to blast falling back to tcoffee. Blast output:\n%s" % (seq.id,r1) )
            mutations, r2 = tcoffee_align( ref, seq )
            if r2 == -1:
                sys.stderr.write( "%s failed to align with tcoffee as well. Tcoffee output:\n%s" % (seq.id,r2) )
        print "%s: Total mutations: %s" % (seq.description, len( mutations ))
        for m in mutations:
            print m
        tcount += 1
