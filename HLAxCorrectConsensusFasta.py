#!/usr/bin/env python
"""

Read fasta file and output a new fasta, where variants have been inserted based on the info in two columns

Author: Rune Moellegaard Friborg <runef@birc.au.dk>
"""

import getopt
import sys
import gzip
from Bio import SeqIO
import vcf
import itertools
import tempfile
import subprocess
import os


##############################################################
####################### Configuration ########################

VERSION="0.09"
UPDATED="2015-11-13"
PID=str(os.getpid())

##############################################################
####################### Classes #################################

class Fasta():
    """
    Reads a fasta files with a single sequence and stores it in self.seq
    """
    def __init__(self, filename):
        sys.stdout.write("Reading '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        if not len(records) == 1:
            sys.stderr.write("Error!: "+str(len(records))+" sequences read. Fasta file '"+ filename +"' must have exactly one sequence!\n")
            sys.exit(1)

        self.val = str(records[0].seq)
        self.id = records[0].id
        sys.stdout.write("done\n")


class Custom():
    """
    Reads a column in a file and stores the list
    """
    def __init__(self, filename, colnum, sep=' ', size=1, has_header=True):
        self.filename = filename
        self.colnum = colnum
        self.data = []

        sys.stdout.write("Parsing '"+filename+"':"+str(colnum)+"...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        if has_header:
            handle.next()

        # Parse lines
        i = 0
        for line in handle:
            try:
                i += 1

                if line == "\n":
                    continue

                l = line.strip().split(sep)
                if size > 1:
                    self.data.append(l[colnum-1:colnum-1+(size)])
                else:
                    self.data.append(l[colnum-1])

            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)
        sys.stdout.write("done\n")
        handle.close()




##############################################################
####################### Main #################################

def main(args):
    # Read Fasta file
    faSequence = Fasta(args.fa_input)

    # Read variants and output new fasta
    posList = Custom(args.var_pos_file, args.var_pos_col).data
    refLenList = Custom(args.var_ref_len_file, args.var_ref_len_col).data
    newSeqList = Custom(args.var_new_seq_file, args.var_new_seq_col).data
    newAltList = Custom(args.var_new_alt_file, args.var_new_alt_col, size=6).data

    vcfInfoList = Custom(args.var_vcf_info_file, 1, size=100).data

    if not (len(posList) == len(refLenList) and len(posList) == len(newSeqList) and len(posList) == len(newAltList) and len(posList) == len(vcfInfoList)):
        sys.stderr.write("Error! Mismatch in length of columns")
        sys.stderr.write("pos: "+str(len(posList))+"\n")
        sys.stderr.write("seq len: "+str(len(refLenList))+"\n")
        sys.stderr.write("new seq: "+str(len(newSeqList))+"\n")
        sys.stderr.write("new alt: "+str(len(newAltList))+"\n")
        sys.stderr.write("vcf info: "+str(len(vcfInfoList))+"\n")
        sys.exit(1)


    # Perform filtering of conflicts
    records = {}
    for i in xrange(len(posList)):
        pos = posList[i]

        if pos == None:
            # entry was removed, because of a conflict
            # skip
            continue

        # type cast to int
        pos = int(pos)

        ref_len = int(refLenList[i])
        for k in range(ref_len):
            if records.has_key(pos+k):
                #possible conflict!
                prev_index = records[pos+k]
                prev_pos = posList[prev_index]
                if prev_pos == None:
                    posList[i] = None
                    continue

                prev_pos = int(prev_pos)
                prev_seq = newSeqList[prev_index]
                cur_seq = newSeqList[i]

                j = pos+k-prev_pos

                seq_len_overlap = min(len(cur_seq)-k, len(prev_seq)-j)

                if j < len(prev_seq) and k < len(cur_seq) and cur_seq[k:k+seq_len_overlap] == prev_seq[j:j+seq_len_overlap]:
                    # Ok
                    pass
                elif (not j < len(prev_seq)) and (not k < len(cur_seq)):
                    sys.stderr.write("WARNING! variant at line:" + str(i+2) + " conflicts with previous variant at line:" + str(prev_index+2) +"\n")
                    sys.stderr.write("  Conflict: deletion VS deletion. Removing both entries.\n")
                    posList[i] = None
                    posList[prev_index] = None
                elif not k < len(cur_seq):
                    sys.stderr.write("WARNING! variant at line:" + str(i+2) + " conflicts with previous variant at line:" + str(prev_index+2) +"\n")
                    sys.stderr.write("  Conflict: '"+prev_seq[j]+"' VS deletion. Removing both entries.\n")
                    posList[i] = None
                    posList[prev_index] = None
                elif not j < len(prev_seq):
                    sys.stderr.write("WARNING! variant at line:" + str(i+2) + " conflicts with previous variant at line:" + str(prev_index+2) +"\n")
                    sys.stderr.write("  Conflict: deletion VS '"+cur_seq[k]+"'. Removing both entries.\n")
                    posList[i] = None
                    posList[prev_index] = None
                else:
                    sys.stderr.write("WARNING! variant at line:" + str(i+2) + " conflicts with previous variant at line:" + str(prev_index+2) +"\n")
                    sys.stderr.write("  Conflict: '"+prev_seq[j:j+seq_len_overlap]+"' VS '"+cur_seq[k:k+seq_len_overlap]+"'. Removing both entries.\n")
                    posList[i] = None
                    posList[prev_index] = None

            else:
                records[pos+k] = i


    sys.stdout.write("Writing '"+args.fa_output+"' and '"+args.vcf_output+"'...")
    ohandle_fasta = open(args.fa_output, "w")
    ohandle_fasta.write(">"+str(faSequence.id)+"\n")
    ohandle_vcf = open(args.vcf_output, "w")

    ohandle_vcf.write('##fileformat=VCFv4.1\n')
    ohandle_vcf.write('##source=HLAxCorrectConsensusFasta.py\n')
    ohandle_vcf.write('##INFO=<ID=PS,Number=1,Type=String,Description="Phasing status">\n')
    ohandle_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    ohandle_vcf.write('#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT '+str(faSequence.id)+'\n')
    #ohandle_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tINFO\n")

    def calc_genotype(method, child_a1, child_a2, var):
        if method == "unassigned":
            return "./."
        elif method == "obvious" or method == "phased":
            if child_a1 == var:
                return "1|1"
            else:
                return "0|0"


    vcf_offset = 0
    region_start = 0
    region_end = 0
    prev_entry = (0, 0, "")
    for entry in zip(posList, refLenList, newSeqList, newAltList, vcfInfoList):
        pos = entry[0]
        if pos == None:
            # Entry removing during filtering. Skip
            continue

        pos = int(pos)
        refLen = int(entry[1])
        newSeq = entry[2]
        newAlt = entry[3]
        vcfInfo = entry[4]

        newAlt = set(newAlt)
        if newSeq in newAlt:
            newAlt.remove(newSeq)

        if pos < region_start:
            prev_pos = int(prev_entry[0])
            prev_refLen = int(prev_entry[1])
            prev_newSeq = prev_entry[2]

            if pos > prev_pos+len(prev_newSeq):
                sys.stderr.write("ERROR! Uncaught conflict in filtering\n")
                sys.stderr.write("  " + str(entry) + " conflicts with previous variant " + str(prev_entry) +"\n")
                sys.exit(1)
            else:
                conflict = False
                # Possible conflict with bases. Check.
                # Example of conflict ('421484', '1', 'AT') conflicts with previous variant ('421480', '5', 'AGTTG')
                # Example of conflict ('50675', '1', 'AG') conflicts with previous variant ('50674', '2', 'T')
                # Example of non-conflict ('5692', '1', 'N') conflicts with previous variant ('5691', '2', 'NN')
                j = 0
                for i in xrange(pos-prev_pos, len(prev_newSeq)):
                    if j < len(newSeq) and prev_newSeq[i] != newSeq[j]:
                        conflict = True
                    j += 1

                if conflict:
                    sys.stderr.write("ERROR! Uncaught conflict in filtering\n")
                    sys.stderr.write("  " + str(entry[4]) + " conflicts with previous variant\n  " + str(prev_entry[4]) +"\n")
                    sys.exit(1)
                else:

                    region_end = prev_pos + len(prev_newSeq)
                    diff = (region_end - pos)

                    # Write VCF entry
                    method = vcfInfo[5]
                    ohandle_vcf.write(str(faSequence.id)+"\t"+str(vcf_offset-diff)+"\told_pos_"+str(pos)+"\t"+str(newSeq)+"\t"+','.join(newAlt)+"\t"+str(vcfInfo[10])+"\tPASS\tPS="+method+"\tGT\t")
                    ohandle_vcf.write(calc_genotype(method, child_a1=vcfInfo[13], child_a2=vcfInfo[14], var=vcfInfo[2]) + "\n")

                    # Write new seq
                    ohandle_fasta.write(newSeq[diff:])
                    vcf_offset += len(newSeq[diff:])

                    # Set region start and skip ref seq
                    region_start = region_end+refLen-diff

        else:
            region_end = pos

            # Write sequence between variants
            ohandle_fasta.write(faSequence.val[region_start:region_end])
            vcf_offset += len(faSequence.val[region_start:region_end])

            # Write VCF entry
            method = vcfInfo[5]
            ohandle_vcf.write(str(faSequence.id)+"\t"+str(vcf_offset)+"\told_pos_"+str(pos)+"\t"+str(newSeq)+"\t"+','.join(newAlt)+"\t"+str(vcfInfo[10])+"\tPASS\tPS="+method+"\tGT\t")
            ohandle_vcf.write(calc_genotype(method, child_a1=vcfInfo[13], child_a2=vcfInfo[14], var=vcfInfo[2]) + "\n")

            # Write new seq
            ohandle_fasta.write(newSeq)
            vcf_offset += len(newSeq)

            # Set region start and skip ref seq
            region_start = region_end+refLen

        # Store entry
        prev_entry = entry

    # Write final region
    ohandle_fasta.write(faSequence.val[region_start:])
    ohandle_fasta.close()
    ohandle_vcf.close()
    sys.stdout.write("done\n")


##############################################################
######################### Help ###############################

def usage():
    print("""HLAxCorrectConsensusFasta version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxCorrectConsensusFasta \
    --fa-input=<fasta file> --var-pos=<file:column> --var-ref-len=<file:column>
    --var-new-seq=<file:column> --var-new-alt=<file:column-offset> --var-vcf-info=<file:column>
    --fa-output=<fasta file> --vcf-output=<vcf file>
""")


class ArgContainer():
    def __init__(self):
        self.fa_input     = ""
        self.fa_output    = ""
        self.vcf_output   = ""
        self.var_pos_file = ""
        self.var_pos_col  = None
        self.var_ref_len_file = ""
        self.var_ref_len_col  = None
        self.var_new_seq_file = ""
        self.var_new_seq_col  = ""
        self.var_new_alt_file = ""
        self.var_new_alt_col  = ""
        self.var_vcf_info_file = ""
        self.var_vcf_info_col  = ""

    def ok(self):
        err = 0
        if not self.fa_input:
            sys.stderr.write("Missing argument: --fa-input\n")
            err = 1
        if not self.var_pos_file:
            sys.stderr.write("Missing file part argument: --var-pos\n")
            err = 1
        if not self.var_pos_col:
            sys.stderr.write("Missing column part argument: --var-pos\n")
            err = 1
        if not self.var_ref_len_file:
            sys.stderr.write("Missing file part argument: --var-ref-len\n")
            err = 1
        if not self.var_ref_len_col:
            sys.stderr.write("Missing column part argument: --var-ref-len\n")
            err = 1
        if not self.var_new_seq_file:
            sys.stderr.write("Missing file part argument: --var-new-seq\n")
            err = 1
        if not self.var_new_seq_col:
            sys.stderr.write("Missing column part argument: --var-new-seq\n")
            err = 1
        if not self.var_new_alt_file:
            sys.stderr.write("Missing file part argument: --var-new-alt\n")
            err = 1
        if not self.var_new_alt_col:
            sys.stderr.write("Missing column part argument: --var-new-alt\n")
            err = 1
        if not self.var_vcf_info_file:
            sys.stderr.write("Missing file part argument: --var-vcf-info\n")
            err = 1
        if not self.var_vcf_info_col:
            sys.stderr.write("Missing column part argument: --var-vcf-info\n")
            err = 1
        if not self.fa_output:
            sys.stderr.write("Missing argument: --fa-output\n")
            err = 1
        if not self.vcf_output:
            sys.stderr.write("Missing argument: --vcf-output\n")
            err = 1
        if err:
            sys.stderr.write("\n")

        return not err



if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "fa-input=", "var-pos=", "var-ref-len=", "var-new-seq=", "var-new-alt=", "var-vcf-info=", "fa-output=", "vcf-output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    args = ArgContainer()
    for o, a in opts:
        if a and a[0] == "=":
            a = a[1:]


        if o == "--fa-input":
            args.fa_input = a
        elif o == "--fa-output":
            args.fa_output = a
        elif o == "--vcf-output":
            args.vcf_output = a
        elif o == "--var-pos":
            try:
                args.var_pos_file, args.var_pos_col = a.split(":")
                args.var_pos_col = int(args.var_pos_col)
            except ValueError:
                pass
        elif o == "--var-ref-len":
            try:
                args.var_ref_len_file, args.var_ref_len_col = a.split(":")
                args.var_ref_len_col = int(args.var_ref_len_col)
            except ValueError:
                pass
        elif o == "--var-new-seq":
            try:
                args.var_new_seq_file, args.var_new_seq_col = a.split(":")
                args.var_new_seq_col = int(args.var_new_seq_col)
            except ValueError:
                pass
        elif o == "--var-new-alt":
            try:
                args.var_new_alt_file, args.var_new_alt_col = a.split(":")
                args.var_new_alt_col = int(args.var_new_alt_col)
            except ValueError:
                pass
        elif o == "--var-vcf-info":
            try:
                args.var_vcf_info_file, args.var_vcf_info_col = a.split(":")
                args.var_vcf_info_col = int(args.var_vcf_info_col)
            except ValueError:
                pass
        elif o == "--help":
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"

    if args.ok():
        try:
            main(args)
        except KeyboardInterrupt:
            sys.exit(1)
    else:
        usage()
        sys.exit()
