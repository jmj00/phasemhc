#!/usr/bin/python
import sys, vcf, gzip, re, time, timeit, pickle, argparse
from Bio import SeqIO
from decimal import Decimal
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

parser = argparse.ArgumentParser(description='''
Find positions of variants in a trio through exact matching of upstream sequences
''')

parser.add_argument('--father_fasta', type = str,
    help="fasta file from father")
parser.add_argument('--mother_fasta', type = str,
    help="fasta file from mother")
parser.add_argument('--child_fasta', type = str,
    help="fasta file from child")
parser.add_argument('--father_vcf', type = str,
    help="vcf file from father")
parser.add_argument('--mother_vcf', type = str,
    help="vcf file from mother")
parser.add_argument('--child_vcf', type = str,
    help="vcf file from child")
parser.add_argument('--father_out', type = str,
    help="output file from father")
parser.add_argument('--mother_out', type = str,
    help="output file from mother")
parser.add_argument('--child_out', type = str,
    help="output file from child")
parser.add_argument('--verbose', action='store_true',
    help="print verbose output")


args = parser.parse_args()


def fileOpener(f):
    if f.endswith('gz'):
        return gzip.open(f)
    else:
        return(open(f))

def make_rc_record(seq):
    """Returns a new SeqRecord sequence with the reverse complement sequence."""
    my_seq = Seq(seq, generic_dna).reverse_complement()
    return SeqRecord(seq = my_seq, \
                 id = "rc", \
                 description = "reverse complement").seq

def getQueryPositions(record_pos, flank_len, record_reference, variant_individual):
    '''Retrieves start and end positions of query sequence from the variant individual'''
    # Handle end cases (for left flanking start and right flanking end)

    if record_pos <= flank_len:
        lqstart = 0
    else:
        lqstart  =   record_pos - (flank_len + 1)

    if record_pos + flank_len + 1 > len(variant_individual):
        rqend = len(variant_individual)
    else:
        rqend = record_pos + len(record_reference) + flank_len + 1

    # Left flanking end and right flanking start will always be defined the same way
    lqend    =   record_pos - 1
    rqstart  =   record_pos + len(record_reference)   # Adjust for length of variant (ref)

    return (lqstart, lqend, rqstart, rqend)

# father_fasta = '/Volumes/GenomeDK/00_HaplotypeAssembly/05_bwbble/datadir/1009-01/1009-01.complete.fa'
# mother_fasta = '/Volumes/GenomeDK/00_HaplotypeAssembly/05_bwbble/datadir/1009-02/1009-02.complete.fa'
# child_fasta = '/Volumes/GenomeDK/00_HaplotypeAssembly/05_bwbble/datadir/1009-06/1009-06.complete.fa'
# vcfdir = '/Volumes/GenomeDK/00_HaplotypeAssembly/05_bwbble/datadir'

father = [seq_record.seq for seq_record in SeqIO.parse(fileOpener(args.father_fasta), 'fasta')][0]
mother = [seq_record.seq for seq_record in SeqIO.parse(fileOpener(args.mother_fasta), 'fasta')][0]
child = [seq_record.seq for seq_record in SeqIO.parse(fileOpener(args.child_fasta), 'fasta')][0]

f_positions = [(record.POS, record.REF, record.ALT) for record in vcf.Reader(open(args.father_vcf, 'r'))]
m_positions = [(record.POS, record.REF, record.ALT) for record in vcf.Reader(open(args.mother_vcf, 'r'))]
c_positions = [(record.POS, record.REF, record.ALT) for record in vcf.Reader(open(args.child_vcf, 'r'))]

flank_len = 40

if args.verbose:
    print >> sys.stderr, "Checking variants from father", len(f_positions)

# fout = '/Volumes/GenomeDK/00_HaplotypeAssembly/05_bwbble/full/1009-01/1009-01.f_rc_out.txt'
# mout = '/Volumes/GenomeDK/00_HaplotypeAssembly/05_bwbble/full/1009-02/1009-02.m_rc_out.txt'
# cout = '/Volumes/GenomeDK/00_HaplotypeAssembly/05_bwbble/full/1009-06/1009-06.c_rc_out.txt'

fout = args.father_out
mout = args.mother_out
cout = args.child_out

with open(fout, 'w') as f:
    f.write('\t'.join(['var.f.pos','ref', 'alt', 'length', 'f.pos', 'm.pos', 'c.pos', '\n']))


for i in range(len(f_positions)):
    if i % 100 == 0:
        print i

    lqstart = getQueryPositions(f_positions[i][0], flank_len, f_positions[i][1], father)[0]
    lqend   = getQueryPositions(f_positions[i][0], flank_len, f_positions[i][1], father)[1]
    lflank = str(father[lqstart:lqend])
    lrcflank = str(make_rc_record(lflank))


    with open(fout, 'a') as f:
        line = "".join(["\t".join([str(f_positions[i][0]), str(f_positions[i][1]), str(f_positions[i][2][0]), str(len(f_positions[i][1]))]), '\t'])
        f.write(line)

    if len([m.end() for m in re.finditer(lflank, str(father))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(father))][0] + 1
        with open(fout, 'a') as f:
            line = "".join([str(currpos),'\t'])
            f.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(father))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(father))][0] - 1) * -1
        with open(fout, 'a') as f:
            line = "".join([str(currpos),'\t'])
            f.write(line)

    else:
        with open(fout, 'a') as f:
            line = "".join(['NA','\t'])
            f.write(line)

    if len([m.end() for m in re.finditer(lflank, str(mother))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(mother))][0] + 1
        with open(fout, 'a') as f:
            line = "".join([str(currpos),'\t'])
            f.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(mother))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(mother))][0] - 1) * -1
        with open(fout, 'a') as f:
            line = "".join([str(currpos),'\t'])
            f.write(line)

    else:
        with open(fout, 'a') as f:
            line = "".join(['NA','\t'])
            f.write(line)

    if len([m.end() for m in re.finditer(lflank, str(child))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(child))][0] + 1
        with open(fout, 'a') as f:
            line = "".join([str(currpos),'\n'])
            f.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(child))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(child))][0] - 1) * -1
        with open(fout, 'a') as f:
            line = "".join([str(currpos),'\n'])
            f.write(line)

    else:
        with open(fout, 'a') as f:
            line = "".join(['NA','\n'])
            f.write(line)


if args.verbose:
    print >> sys.stderr, "Checking variants from mother", len(m_positions)

with open(mout, 'w') as m:
    m.write('\t'.join(['var.m.pos','ref', 'alt', 'length', 'f.pos', 'm.pos', 'c.pos', '\n']))

for i in range(len(m_positions)):
    if i % 100 == 0:
        print i

    lqstart = getQueryPositions(m_positions[i][0], flank_len, m_positions[i][1], mother)[0]
    lqend   = getQueryPositions(m_positions[i][0], flank_len, m_positions[i][1], mother)[1]
    lflank = str(mother[lqstart:lqend])

    with open(mout, 'a') as m:
        line = "".join(["\t".join([str(m_positions[i][0]), str(m_positions[i][1]), str(m_positions[i][2][0]), str(len(m_positions[i][1]))]), '\t'])
        m.write(line)

    if len([m.end() for m in re.finditer(lflank, str(father))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(father))][0] + 1
        with open(mout, 'a') as m:
            line = "".join([str(currpos),'\t'])
            m.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(father))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(father))][0] - 1) * -1
        with open(mout, 'a') as m:
            line = "".join([str(currpos),'\t'])
            m.write(line)

    else:
        with open(mout, 'a') as m:
            line = "".join(['NA','\t'])
            m.write(line)

    if len([m.end() for m in re.finditer(lflank, str(mother))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(mother))][0] + 1
        with open(mout, 'a') as m:
            line = "".join([str(currpos),'\t'])
            m.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(mother))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(mother))][0] - 1) * -1
        with open(mout, 'a') as m:
            line = "".join([str(currpos),'\t'])
            m.write(line)

    else:
        with open(mout, 'a') as m:
            line = "".join(['NA','\t'])
            m.write(line)

    if len([m.end() for m in re.finditer(lflank, str(child))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(child))][0] + 1
        with open(mout, 'a') as m:
            line = "".join([str(currpos),'\n'])
            m.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(child))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(child))][0] - 1) * -1
        with open(mout, 'a') as m:
            line = "".join([str(currpos),'\n'])
            m.write(line)

    else:
        with open(mout, 'a') as m:
            line = "".join(['NA','\n'])
            m.write(line)

if args.verbose:
    print >> sys.stderr, "Checking variants from child", len(c_positions)
with open(cout, 'w') as c:
    c.write('\t'.join(['var.c.pos','ref', 'alt', 'length', 'f.pos', 'm.pos', 'c.pos', '\n']))

for i in range(len(c_positions)):
    if i % 100 == 0:
        print i

    lqstart = getQueryPositions(c_positions[i][0], flank_len, c_positions[i][1], child)[0]
    lqend   = getQueryPositions(c_positions[i][0], flank_len, c_positions[i][1], child)[1]
    lflank = str(child[lqstart:lqend])

    with open(cout, 'a') as c:
        line = "".join(["\t".join([str(c_positions[i][0]), str(c_positions[i][1]), str(c_positions[i][2][0]), str(len(c_positions[i][1]))]), '\t'])
        c.write(line)


    if len([m.end() for m in re.finditer(lflank, str(father))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(father))][0] + 1
        with open(cout, 'a') as c:
            line = "".join([str(currpos),'\t'])
            c.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(father))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(father))][0] - 1) * -1
        with open(cout, 'a') as c:
            line = "".join([str(currpos),'\t'])
            c.write(line)

    else:
        with open(cout, 'a') as c:
            line = "".join(['NA','\t'])
            c.write(line)

    if len([m.end() for m in re.finditer(lflank, str(mother))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(mother))][0] + 1
        with open(cout, 'a') as c:
            line = "".join([str(currpos),'\t'])
            c.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(mother))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(mother))][0] - 1) * -1
        with open(cout, 'a') as c:
            line = "".join([str(currpos),'\t'])
            c.write(line)

    else:
        with open(cout, 'a') as c:
            line = "".join(['NA','\t'])
            c.write(line)

    if len([m.end() for m in re.finditer(lflank, str(child))]) == 1:
        currpos = [m.end() for m in re.finditer(lflank, str(child))][0] + 1
        with open(cout, 'a') as c:
            line = "".join([str(currpos),'\n'])
            c.write(line)

    elif len([m.end() for m in re.finditer(lrcflank, str(child))]) == 1:
        currpos = ([m.end() for m in re.finditer(lrcflank, str(child))][0] - 1) * -1
        with open(cout, 'a') as c:
            line = "".join([str(currpos),'\n'])
            c.write(line)

    else:
        with open(cout, 'a') as c:
            line = "".join(['NA','\n'])
            c.write(line)


with open(fout, 'a') as f:
    f.write('\n')
with open(mout, 'a') as m:
    m.write('\n')
with open(cout, 'a') as c:
    c.write('\n')





