#!/usr/bin/python
import re, sys, collections, numpy, gzip, vcf, argparse, subprocess, time
from subprocess import call
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Concatenate scaffold fastas in correct order, by using alignmentblocks from MAF file. Adjust VCFs accordingly')

parser.add_argument('-m', '--maf', action='store', dest='maf', type = str,
                    help='maf file of alignment of scaffolds to reference')

parser.add_argument('-s', '--scaffold-list', action='store', dest='scaffold_list', type = str,
                    help='list of scaffolds')

parser.add_argument('-i', '--in-prefix', action='store', dest='iprefix', type = str,
                    help='input prefix of fasta and vcf file')

parser.add_argument('-o', '--out-prefix', action='store', dest='oprefix', type = str,
                    help='output prefix of fasta and vcf file')

parser.add_argument('-c', '--complete-prefix', action='store', dest='cprefix', type = str,
                    help='prefix of comcplete (concatenated) fasta and vcf file')

parser.add_argument('--HLA-F', action='store', dest='HLAF_fasta', type = str,
                    help='HLA-F fasta file')

parser.add_argument('--KIFC1', action='store', dest='KIFC1_fasta', type = str,
                    help='KIFC1 fasta file')

parser.add_argument('--blast-dir', action='store', dest='blast_dir', type = str,
                    help='Directory of BLAST binaries')


args = parser.parse_args()


def makeblastdb(fasta, db_type, db_name, db_title, blast_dir):
    """Build blast database from local fasta file"""
    blast_executable = "".join([str(blast_dir), "makeblastdb"])

    f = gzip.open(fasta)
    format_cmd = '%s -in - -dbtype %s -out %s -title %s' % (blast_executable, db_type, db_name, db_title)
    process = subprocess.Popen(format_cmd,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    shell=True)
    process.communicate(f.read())


def blastn(db_name, query_fasta, out_file, blast_dir):
    """Perform nucleotide blast (blastn) using local database"""
    print "Blasting ", query_fasta, "against ", db_name, "\n"
    blast_executable = "".join([str(blast_dir), "blastn"])

    format_cmd = '%s -db %s -query %s -out %s' % (blast_executable, db_name, query_fasta, out_file)
    process = subprocess.Popen(format_cmd,
                                    shell=True)
    print "Sleeping for 40 seconds. Issues encountered: ", time.sleep(40)


def parseBlast(blast_out):
    """Parses output of BLAST to get start position of best hit"""
    with open(blast_out) as hits:
        query_start = []
        sbjct_start = []
        first_match = 1
        for l in hits:
            if l.startswith("Query"):
                if not first_match == 1:
                    query_start.append(l.split()[1])
                first_match = 0
            if l.startswith("Sbjct"):
                sbjct_start.append(l.split()[1])

        if sbjct_start and query_start:
            return int(sbjct_start[0]) - int(query_start[0]) + 1
        else:
            print "No hits"
def median(lst):
    return numpy.median(numpy.array(lst))

def fileOpener(f):
    if f.endswith('gz'):
        return gzip.open(f)
    else:
        return(open(f))

def mafReader(maf, scaffold_list):
    """Given a maf file and list of scaffolds, parse the maf records and create dictionary with alignment info"""
    records = {}
    scaffolds = getScaffoldList(scaffold_list)
    x = 0
    with fileOpener(maf) as f:
        switch = 0
        for l in f:

            if re.search('chr6_', l):
                switch = 1
                chr_line = l.strip().split()

            if re.search('scaffold', l) and switch == 1:
                entry = l.strip().split()
                scaffold = entry[1].split('_')[1]
                switch = 0

                if scaffold not in records and scaffold in scaffolds:
                    records[scaffold] = {'sense': [], 'antisense': [], 'seq': [], 'chr_block': [], 'median': 0, 'reverse_complement': 0}
                    records[scaffold]['seq'] = entry[6]

                if scaffold in scaffolds:
                    records[scaffold]['chr_block'].append(int(chr_line[2]))

                if entry[4] == '+' and scaffold in scaffolds:
                    records[scaffold]['sense'].append(entry[3])

                elif entry[4] == '-' and scaffold in scaffolds:
                    records[scaffold]['antisense'].append(entry[3])


        for scaffold in records:

            sense = 0
            antisense = 0
            records[scaffold]['median'] = median(records[scaffold]['chr_block'])

            sense_antisense = records[scaffold]['sense'] + records[scaffold]['antisense']
            max_block = max(map(int, sense_antisense))
            print "Found largest maximum alignment block to be :", int(max_block)
            print "Filtering out small alignment blocks"
            sense_list = [b for b in records[scaffold]['sense'] if int(b) >= max_block/2]
            antisense_list = [b for b in records[scaffold]['antisense'] if int(b) >= max_block/2]
            sense = sum(map(int, sense_list))
            antisense = sum(map(int, antisense_list))

            if sense >= antisense:
                records[scaffold]['reverse_complement'] = 0

            else:
                records[scaffold]['reverse_complement'] = 1

        return records

def getScaffoldList(scaffold_list):
    """Read in list of scaffolds"""
    with open(scaffold_list) as f:
        scaffolds = [line.strip() for line in f]

    return scaffolds

def getOrder(dictionary):
    """Get the sorted order of scaffolds by median of alignment start positions"""
    M = []
    for scaffold in dictionary:
        M.append([scaffold, dictionary[scaffold]['median']])

    M.sort(key = lambda row: row[1])

    return [x[0] for x in M]


def reverseComplement(record):
    """Takes a sequence (str) as input and returns a new SeqRecord with the reverse complement sequence."""
    my_seq = Seq(record, generic_dna).reverse_complement()

    return my_seq

def makeSeqRecord(string, record_name, record_id, record_description):
    """Create SeqRecord object"""
    return SeqRecord(name = record_name, \
                 seq = Seq(string), \
                 id = record_id, \
                 description = record_description)

def adjustRecords(dictionary, order, in_prefix, out_prefix, complete_prefix, HLAF_fasta, KIFC1_fasta, blast_dir):
    """Adjust coordinates of vcf file and reverse complement sequence for reverse complemented scaffolds in vcf and fasta, BLAST and trim"""

    for i, scaffold in enumerate(order):

        ifasta = ".".join([in_prefix, str(scaffold), 'fa.gz'])
        ofasta = ".".join([out_prefix, str(scaffold), 'fa.gz'])
        complete_fasta = ".".join([complete_prefix, 'fa'])
        trio_id = complete_prefix.split("/")[9].split("-")[0]

        sequence = str([seq_record.seq for seq_record in SeqIO.parse(fileOpener(ifasta), 'fasta')][0])
        scaffold_length = len(sequence)

        if i == 0 and i == len(order) - 1:
            # If only one scaffolds is found
            if dictionary[scaffold]['reverse_complement'] == 1:
                rc_fasta = "".join([ifasta.split('fa.gz')[0], "rc.fa.gz"])
                seq_rc = reverseComplement(sequence)
                with gzip.open(rc_fasta, 'w') as out:
                    seq_rc_out = makeSeqRecord(str(seq_rc), scaffold+"_rc", scaffold, "reverse_complemented")
                    out.write(seq_rc_out.format('fasta'))

                db_name = ".".join([in_prefix, str(scaffold), "rc.fa"])
                db_type = 'nucl'
                db_title = str(scaffold)
                query_fasta = HLAF_fasta
                HLAF_out = ".".join([in_prefix, str(scaffold), "HLA-F", "out"])
                makeblastdb(rc_fasta, db_type, db_name, db_title, blast_dir)
                blastn(db_name, query_fasta, HLAF_out, blast_dir)
                start = parseBlast(HLAF_out)
                query_fasta = KIFC1_fasta
                KIFC1_out = ".".join([in_prefix, str(scaffold), "KIFC1", "out"])
                blastn(db_name, query_fasta, KIFC1_out, blast_dir)
                end = parseBlast(KIFC1_out)
                edge = end + 18387 + 1000 # (+ length of KIFC1 gene + 1kb)

                if start <= 1:
                    offset = 1
                else:
                    offset = -(start - 1000000 - 1)

                trimmed_seq = makeSeqRecord(str(seq_rc[int(-offset-1):edge]), scaffold+"_trimmed", scaffold, "trimmed")

            else:
                db_name = ".".join([in_prefix, str(scaffold), "rc.fa"])
                db_type = 'nucl'
                db_title = str(scaffold)
                query_fasta = HLAF_fasta
                HLAF_out = ".".join([in_prefix, str(scaffold), "HLA-F", "out"])
                makeblastdb(ifasta, db_type, db_name, db_title, blast_dir)
                blastn(db_name, query_fasta, HLAF_out, blast_dir)
                start = parseBlast(HLAF_out)
                query_fasta = KIFC1_fasta
                KIFC1_out = ".".join([in_prefix, str(scaffold), "KIFC1", "out"])
                blastn(db_name, query_fasta, KIFC1_out, blast_dir)
                end = parseBlast(KIFC1_out)
                edge = end + 18387 + 1000 # (+ length of KIFC1 gene + 1kb)

                if start <= 1:
                    offset = 1
                else:
                    offset = -(start - 1000000 - 1)

                trimmed_seq = makeSeqRecord(str(sequence[int(-offset-1):edge]), scaffold+"_trimmed", scaffold, "trimmed")

            with gzip.open(ofasta, 'w') as out:
                out.write(trimmed_seq.format('fasta'))

            final_sequence = str(trimmed_seq.seq)

            with open(complete_fasta, 'w') as out:
                final_seq = makeSeqRecord(str(final_sequence), "complete", trio_id, "complete hla")

                out.write(final_seq.format('fasta'))

        elif i == 0:

            if dictionary[scaffold]['reverse_complement'] == 1:

                rc_fasta = "".join([ifasta.split('fa.gz')[0], "rc.fa.gz"])
                seq_rc = reverseComplement(sequence)
                with gzip.open(rc_fasta, 'w') as out:
                    seq_rc_out = makeSeqRecord(str(seq_rc), scaffold+"_rc", scaffold, "reverse_complemented")
                    out.write(seq_rc_out.format('fasta'))

                db_name = ".".join([in_prefix, str(scaffold), "rc.fa"])
                db_type = 'nucl'
                db_title = str(scaffold)
                query_fasta = HLAF_fasta
                HLAF_out = ".".join([in_prefix, str(scaffold), "HLA-F", "out"])
                makeblastdb(rc_fasta, db_type, db_name, db_title, blast_dir)
                blastn(db_name, query_fasta, HLAF_out, blast_dir)
                start = parseBlast(HLAF_out)

                if start <= 1:
                    offset = 1
                else:
                    offset = -(start - 1000000 - 1)

                trimmed_seq = makeSeqRecord(str(seq_rc[int(-offset-1):int(scaffold_length)])+'N', scaffold+"_trimmed", scaffold, "trimmed")

            else:
                db_name = ".".join([in_prefix, str(scaffold), "fa"])
                db_type = 'nucl'
                db_title = str(scaffold)
                query_fasta = HLAF_fasta
                HLAF_out = ".".join([in_prefix, str(scaffold), "HLA-F", "out"])
                makeblastdb(ifasta, db_type, db_name, db_title, blast_dir)
                blastn(db_name, query_fasta, HLAF_out, blast_dir)
                start = parseBlast(HLAF_out)

                if start <= 1:
                    offset = 1
                else:
                    offset = -(start - 1000000 - 1)
                trimmed_seq = makeSeqRecord(str(sequence[int(-offset-1):int(scaffold_length)])+'N', scaffold+"_trimmed", scaffold, "trimmed")

            with gzip.open(ofasta, 'w') as out:
                out.write(trimmed_seq.format('fasta'))

            final_sequence = str(trimmed_seq.seq)

        elif i == len(order)-1:

            trimmed_seq = makeSeqRecord(str(sequence), scaffold+"_trimmed", scaffold, "trimmed")
            if dictionary[scaffold]['reverse_complement'] == 1:

                rc_fasta = ".".join([ifasta.split('.fa.gz')[0], "rc.fa.gz"])
                seq_rc = reverseComplement(sequence)
                with gzip.open(rc_fasta, 'w') as out:
                    seq_rc_out = makeSeqRecord(str(seq_rc), scaffold+"_rc", scaffold, "reverse_complemented")
                    out.write(seq_rc_out.format('fasta'))

                db_name = ".".join([in_prefix, str(scaffold), "rc.fa"])
                db_type = 'nucl'
                db_title = str(scaffold)
                query_fasta = KIFC1_fasta
                KIFC1_out = ".".join([in_prefix, str(scaffold), "KIFC1", "out"])
                makeblastdb(rc_fasta, db_type, db_name, db_title, blast_dir)
                blastn(db_name, query_fasta, KIFC1_out, blast_dir)
                end = parseBlast(KIFC1_out)
                edge = end + 18387 + 1000 # (+ length of KIFC1 gene + 1kb)
                trimmed_seq = makeSeqRecord(str(seq_rc[:edge]), scaffold+"_trimmed", scaffold, "trimmed")

            else:
                db_name = ".".join([in_prefix, str(scaffold), "fa"])
                db_type = 'nucl'
                db_title = str(scaffold)
                query_fasta = KIFC1_fasta
                KIFC1_out = ".".join([in_prefix, str(scaffold), "KIFC1", "out"])
                makeblastdb(ifasta, db_type, db_name, db_title, blast_dir)
                blastn(db_name, query_fasta, KIFC1_out, blast_dir)
                end = parseBlast(KIFC1_out)

                edge = end + 18387 + 1000 # (+ length of KIFC1 gene + 1kb)
                trimmed_seq = makeSeqRecord(str(sequence[:edge]), scaffold+"_trimmed", scaffold, "trimmed")

            with gzip.open(ofasta, 'w') as out:
                out.write(trimmed_seq.format('fasta'))

            final_sequence = final_sequence + str(trimmed_seq.seq)

            with open(complete_fasta, 'w') as out:
                final_seq = makeSeqRecord(str(final_sequence), "complete", trio_id, "complete hla")

                out.write(final_seq.format('fasta'))

        else:
            if dictionary[scaffold]['reverse_complement'] == 1:
                final_sequence = final_sequence + str(reverseComplement(sequence)) + 'N'
            else:
                final_sequence = final_sequence + sequence + 'N'


        ivcf = ".".join([in_prefix, str(scaffold), 'vcf'])
        ovcf = ".".join([out_prefix, str(scaffold), 'vcf'])

        vcf_reader = vcf.Reader(open(ivcf, 'r'))
        vcf_writer = vcf.Writer(open(ovcf, 'w'), vcf_reader)

        if dictionary[scaffold]['reverse_complement'] == 1:

            for record in vcf_reader:
                record.CHROM = trio_id
                record.ID = record.ID + "_" + scaffold
                record.POS = record.POS - 2
                if i == 0:
                    record.REF = str(reverseComplement(record.REF))
                    record.ALT = [str(reverseComplement(str(record.ALT[0])))]
                    record.POS = scaffold_length - record.POS - (len(record.REF)-1)
                    if record.POS <= scaffold_length + offset:
                        record.POS = record.POS-abs(offset)
                        if record.POS > 0:
                            vcf_writer.write_record(record)

                elif i == len(order) - 1:

                    record.REF = str(reverseComplement(record.REF))
                    record.ALT = [str(reverseComplement(str(record.ALT[0])))]
                    record.POS = offset - record.POS - (len(record.REF)-1)
                    record.POS = record.POS + scaffold_length

                    if record.POS > 0 and record.POS <= offset + edge:
                        vcf_writer.write_record(record)

                else:
                    record.REF = str(reverseComplement(record.REF))
                    record.ALT = [str(reverseComplement(str(record.ALT[0])))]
                    record.POS = offset - record.POS - (len(record.REF)-1)
                    record.POS = record.POS + scaffold_length

                    if record.POS > 0:
                        vcf_writer.write_record(record)
        else:
            if i == len(order) - 1:
                for record in vcf_reader:
                    record.CHROM = trio_id
                    record.ID = record.ID + "_" + scaffold
                    record.POS = record.POS + offset + 1
                    if record.POS >= 0 and record.POS <= offset + edge:
                        vcf_writer.write_record(record)

            else:
                for record in vcf_reader:
                    record.CHROM = trio_id
                    record.ID = record.ID + "_" + scaffold
                    record.POS = record.POS + offset + 1
                    if record.POS >= 0:
                        vcf_writer.write_record(record)

        offset += scaffold_length + 1



def concatenateVCFs(file_list, complete_vcf):
    """Concatenates a list of vcf files"""
    with open(file_list[0]) as first:

        first_line = 1
        for l in first:
            if re.search("#", l):
                if first_line == 1:
                    with open(complete_vcf, 'w') as cout:
                        cout.write(l)
                        first_line = 0
                else:
                    with open(complete_vcf, 'a') as cout:
                        cout.write(l)
            else:
                with open(complete_vcf, 'a') as cout:
                    cout.write(l)

    for r in file_list[1:]:
        with open(r) as remaining:
            for l in remaining:
                if not l.startswith("#"):
                    with open(complete_vcf, 'a') as cout:
                        cout.write(l)


MAF = mafReader(args.maf, args.scaffold_list)
order = getOrder(MAF)
print "______________________________________________________"
print "Order of scaffolds is :", [order[s] for s in range(len(order))]
print "-----------------------"
print "Scaffolds that were reverse complemented:", [s for s in MAF if MAF[s]['reverse_complement'] == 1]
print "Scaffolds that were not reverse complemented:", [s for s in MAF if MAF[s]['reverse_complement'] == 0]
print "______________________________________________________"

adjustRecords(MAF, order, args.iprefix, args.oprefix, args.cprefix, args.HLAF_fasta, args.KIFC1_fasta, args.blast_dir)
vcfList = [".".join([args.oprefix, scaffold, "vcf"]) for scaffold in order]
cvcf = ".".join([args.cprefix, 'vcf'])
concatenateVCFs(vcfList, cvcf)






