#!/usr/bin/env python
"""

Read input of format:
var.f.pos	ref	alt	length	f.pos   m.pos   c.pos
60		T	C	1	60      NA      NA
2305		G	A	1       2305    NA      NA
5181		AAAAAAA	AAAAAA	7    	5181    1541    1716
6817		T	C	1       6817    3170    3343

Read the matching VCF and Fasta for the data above.

For every entry in the primary format, find the variant in either the VCF file or the fasta.

Output final results for phasing.

Author: Rune Moellegaard Friborg <runef@birc.au.dk>
"""


import getopt
import sys
import gzip
from Bio import SeqIO
import vcf

STATS = {
    'OUTPUT_VARIANTS_TOTAL_USED':0,
    'OUTPUT_VARIANTS_MULTIALLELIC_SKIPPED':0,
    'OUTPUT_VARIANTS_NEGATIVE_SKIPPED':0,
    'OUTPUT_VARIANTS_CHILD_NA_SKIPPED':0
}

def print_stats():
    global STATS
    keys = STATS.keys()
    keys.sort()
    for k in keys:
        print('{:<40}:{:>10}'.format(k, STATS[k]))



##############################################################
####################### Configuration ########################

VERSION="0.05"
UPDATED="2015-06-16"

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
        
        self.val = records[0].seq
        sys.stdout.write("done\n")


class VCF():
    """
    Reads a VCF file and stores all records by position in self.dict
    """
    def __init__(self, filename):
        sys.stdout.write("Reading '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")
        vcf_reader = vcf.Reader(handle)
        vcfDict = {}
        for record in vcf_reader:
            vcfDict[record.POS] = record
        handle.close()

        self.val = vcfDict
        sys.stdout.write("done\n")


class JMJProcess():
    """
    Descripe the process.
    """
    def __init__(self, vcf, seq):
        self.child_vcf = vcf.val
        self.child_seq = seq.val
        self.custom_variants = []

    def preprocess(self, filename):
        """
        For every variant in the mother and father identify the variants in the child
        
        Use this method for the father and mother
        """
        new_variants= {}

        sys.stdout.write("Preprocessing '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        # Skip header
        handle.next()

        # Parse lines
        i = 0
        for line in handle:
            try:
                i += 1

                if line == "\n":
                    continue

                l = line.strip().split('\t')
                pos, ref, alt, length, f_pos, m_pos, c_pos = l
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)            

            # Skip if NA
            if (c_pos == 'NA'):
                STATS['OUTPUT_VARIANTS_CHILD_NA_SKIPPED'] += 1
                continue

            # Type cast to int
            pos = int(pos)    
            c_pos = int(c_pos)

            if (pos < 0 or c_pos < 0):
                STATS['OUTPUT_VARIANTS_NEGATIVE_SKIPPED'] += 1
                continue

            # Extract variants
            v = [ref, alt]

            if self.child_vcf.has_key(c_pos):
                # Child has a variant on the same position in it's VCF                
                child_record = self.child_vcf[c_pos]

                cv = [child_record.REF]
                cv.extend([x.sequence for x in child_record.ALT])

                if new_variants.has_key(c_pos):
                    # Variant conflict, Ignore variant!
                    new_variants[c_pos] = None
                else:
                    new_variants[c_pos] = ((f_pos, m_pos), cv)
            else:
                # Search fasta

                # Correct for one-indexing TODO!!! WARNING. Must be changed when vcf files have been corrected
                # fasta_c_pos = c_pos-1
                fasta_c_pos = c_pos

                # Get variants
                cv = []
                for item in v:
                    cv.append(str(self.child_seq[fasta_c_pos:(fasta_c_pos+len(item))]))

                if new_variants.has_key(c_pos):
                    # Variant conflict, Ignore variant!
                    new_variants[c_pos] = None
                else:
                    new_variants[c_pos] = ((f_pos, m_pos), cv)
        
        # Add newly found custom variants
        self.custom_variants.append(new_variants)
        
        sys.stdout.write("done\n")
        handle.close()
        
    def process(self, filename, ofilename, parent_vcf_obj, parent_seq_obj):

        parent_vcf = [x.val for x in parent_vcf_obj]
        parent_seq = [x.val for x in parent_seq_obj]

        sys.stdout.write("Processing '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        # Open output file
        if ofilename[-2:] == 'gz':
            ohandle = gzip.open(ofilename, "w")
            ohandle_skipped = gzip.open(ofilename + ".skipped", "w")
        else:
            ohandle = open(ofilename, "w")
            ohandle_skipped = open(ofilename + ".skipped", "w")


        # Skip header
        handle.next()

        
        # Parse lines
        i = 0
        data = {}
        for line in handle:
            try:
                i += 1

                if line == "\n":
                    continue

                l = line.strip().split('\t')
                pos, ref, alt, length, f_pos, m_pos, c_pos = l
                
                data[int(pos)] = (ref, alt, length, f_pos, m_pos)
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1) 

        # Find max position
        max_pos = 0
        for seq in parent_seq + [self.child_seq]:
            if len(seq) > max_pos:
                max_pos = len(seq)

        # Get Chromoson id
        PRE_CHROM = self.child_vcf.values()[0].CHROM[:4] 
        CHILD_CHROM = self.child_vcf.values()[0].CHROM

        # Write header
        ohandle.write('##fileformat=VCFv4.1\n')
        ohandle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        ohandle.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Position">\n')
        #ohandle.write('##INFO=<ID=VT,Number=1,Type=String,Description="indicates what type of variant the line represents">')
        ohandle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}-01\t{0}-02\t{1}\n".format(PRE_CHROM, CHILD_CHROM))

        
        # Do a full sweep of all positions and construct the resulting record
        for i in xrange(0, max_pos):
            
            CHILD_POS = None
            CHILD_REF = None
            CHILD_ALT = None
            PARENT_POS = [None, None]
            PARENT_REF = [None, None]
            PARENT_ALT = [None, None]
                        
            if data.has_key(i):
                ref, alt, length, f_pos, m_pos = data[i]

                CHILD_POS = i
                CHILD_REF = ref
                CHILD_ALT = alt
                PARENT_POS[0] = f_pos
                PARENT_POS[1] = m_pos

                #sys.stdout.write("X:"+str(CHILD_POS)+":"+ref+","+alt+":f_pos="+f_pos+",m_pos="+m_pos+"\n")

                for j,p,v,s in zip([0,1], PARENT_POS, parent_vcf, parent_seq):
                    
                    if p == 'NA':
                        continue

                    if v.has_key(int(p)):
                        parent_record = v[int(p)]
                        
                        PARENT_REF[j] = parent_record.REF
                        PARENT_ALT[j] = parent_record.ALT[0].sequence

                        if len(parent_record.ALT) != 1:
                            sys.stderr.write("Error! Parent VCF record with multiple ALT variants is not allowed!\n")
                            sys.exit(1)

                        #sys.stdout.write("  parent_vcf:"+str(j)+":"+str([PARENT_REF[j], PARENT_ALT[j]])+"\n")

                    else:
                        # entry not found in VCF

                        # Search fasta

                        # Correct for one-indexing TODO!!! WARNING. Must be changed when vcf files have been corrected
                        # fasta_c_pos = c_pos-1
                        fasta_p = int(p)

                        PARENT_REF[j] = str(s[fasta_p:(fasta_p+len(CHILD_REF))])
                        PARENT_ALT[j] = str(s[fasta_p:(fasta_p+len(CHILD_ALT))])

                        #sys.stdout.write("  parent_fasta:"+str(j)+":"+str([PARENT_REF[j], PARENT_ALT[j]])+"\n")

            p = 0
            for D in self.custom_variants:
                if D.has_key(i) and D[i] != None:
                    #sys.stdout.write("Y:"+str(i)+":p="+str(p)+":"+str(D[i])+"\n")

                    CHILD_POS = i

                    # Always select the shortest variant from two variant from same position found in both parents
                    if CHILD_REF == None or len(D[i][1][0]) < len(CHILD_REF):
                        CHILD_REF = D[i][1][0]                        
                    if CHILD_ALT == None or len(D[i][1][1]) < len(CHILD_ALT):
                        CHILD_ALT = D[i][1][1]

                    # Fetch variants from both parents
                    for sub_p in range(2):
                        p_pos = D[i][0][sub_p]
                        
                        # Skip if NA
                        if p_pos == 'NA':
                            continue

                        # Typecast to int
                        p_pos = int(p_pos)

                        # Fetch from VCF
                        if parent_vcf[sub_p].has_key(p_pos):
                            parent_record = parent_vcf[sub_p][p_pos]

                            PARENT_POS[sub_p] = p_pos
                            PARENT_REF[sub_p] = parent_record.REF
                            PARENT_ALT[sub_p] = parent_record.ALT[0].sequence

                            if len(parent_record.ALT) != 1:
                                sys.stderr.write("Error! Parent VCF record with multiple ALT variants is not allowed!\n")
                                sys.exit(1)

                            #sys.stdout.write("  parent_vcf:"+str(sub_p)+":"+str([PARENT_REF[sub_p], PARENT_ALT[sub_p]])+"\n")
                        elif p != sub_p:

                            # entry not found in VCF
                            # Search fasta

                            # Skip searching fasta if already set
                            if PARENT_POS[sub_p] != None:
                                continue

                            s = parent_seq[sub_p]

                            # Correct for one-indexing TODO!!! WARNING. Must be changed when vcf files have been corrected
                            # fasta_c_pos = c_pos-1
                            fasta_pos = int(p_pos) 

                            PARENT_POS[sub_p] = p_pos
                            PARENT_REF[sub_p] = str(s[fasta_pos:(fasta_pos+len(CHILD_REF))])
                            PARENT_ALT[sub_p] = str(s[fasta_pos:(fasta_pos+len(CHILD_ALT))])

                            #sys.stdout.write("  parent_fasta:"+str(sub_p)+":"+str([PARENT_REF[sub_p], PARENT_ALT[sub_p]])+"\n")

                        else:
                            # If value was not found in VCF, then p_pos is marked invalid and can be skipped
                            # The reason being that the indexes are generated from the VCF file.
                            pass
                        
                p += 1

            if CHILD_POS != None:

                # CHILD_POS, CHILD_REF, CHILD_ALT
                # PARENT_POS[0,1], PARENT_REF[0,1], PARENT_ALT[0,1]

                if ((CHILD_POS != None and CHILD_POS < 0) or
                    (PARENT_POS[0] != None and PARENT_POS[0] < 0) or
                    (PARENT_POS[1] != None and PARENT_POS[1] < 0)):
                    STATS['OUTPUT_VARIANTS_NEGATIVE_SKIPPED'] += 1
                    continue
                
                # Construct variant list
                V = [CHILD_REF, CHILD_ALT]
                for p in range(2):
                    if PARENT_REF[p] and PARENT_ALT[p]:                        
                        V.append(PARENT_REF[p])
                        V.append(PARENT_ALT[p])
                
                def uniq(input):
                    output = []
                    for x in input:
                        if x not in output:
                            output.append(x)
                    return output

                # Make variant list into a unique list
                V = uniq(V)

                # Skip entry if the variant list has more or less than two variants.
                if (len(V) == 2):
                    h = ohandle
                    STATS['OUTPUT_VARIANTS_TOTAL_USED']+= 1
                else:
                    h = ohandle_skipped
                    STATS['OUTPUT_VARIANTS_MULTIALLELIC_SKIPPED']+= 1

                h.write("{0}\t{1}\tvariant_{1}\t{2}\t{3}\t255\tPASS\t\tGT:PS".format(PRE_CHROM, CHILD_POS, V[0], ",".join(V[1:])))
                
                # Output genotype for father and mother
                for p in range(2):
                    if PARENT_REF[p] and PARENT_ALT[p]: 
                        h.write("\t{0}/{1}:{2}".format(V.index(PARENT_REF[p]), V.index(PARENT_ALT[p]), PARENT_POS[p]))
                    else:
                        h.write("\t./.:0")

                # Output genotype for child
                h.write("\t{0}/{1}:{2}\n".format(V.index(CHILD_REF), V.index(CHILD_ALT), CHILD_POS))


        sys.stdout.write("done\n")
        ohandle.close()
        ohandle_skipped.close()
        handle.close()
    

        
        
##############################################################
####################### Main #################################

def main(args):

    # Read vcf file
    vcfRecords = []
    for filename in [args.f_vcf, args.m_vcf, args.c_vcf]:
        vcfRecords.append(VCF(filename))

    # Read Fasta files
    faSequence = []
    for filename in [args.f_fa, args.m_fa, args.c_fa]:
        faSequence.append(Fasta(filename))

    # The main object for storing and computing the phasing values.
    jmj = JMJProcess(vcfRecords[2], faSequence[2])

    # Preprocessing the data for father and mother
    for filename in [args.f_input, args.m_input]:
        jmj.preprocess(filename)
    
    #Processing the child
    jmj.process(args.c_input, args.c_output, vcfRecords[0:2], faSequence[0:2])


    print_stats()
    



##############################################################
######################### Help ###############################

def usage():
    print("""HLAxPreparePhasing version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxPreparePhasing --f-vcf=<file> --f-fa=<fasta file> --f-input=<input file> \
    --m-vcf=<file> --m-fa=<fasta file> --m-input=<input file> \
    --c-vcf=<file> --c-fa=<fasta file> --c-input=<input file> --c-output=<output file>
""")


class ArgContainer():
    def __init__(self):
        self.f_vcf    = ""
        self.f_fa     = ""
        self.f_input  = ""
        self.m_vcf    = ""
        self.m_fa     = ""
        self.m_input  = ""
        self.c_vcf    = ""
        self.c_fa     = ""
        self.c_input  = ""
        self.c_output  = ""

    def ok(self):
        err = 0
        if not self.f_vcf:
            sys.stderr.write("Missing argument: --f-vcf\n")
            err = 1
        if not self.f_fa:
            sys.stderr.write("Missing argument: --f-fa\n")
            err = 1
        if not self.f_input:
            sys.stderr.write("Missing argument: --f-input\n")
            err = 1
        if not self.m_vcf:
            sys.stderr.write("Missing argument: --m-vcf\n")
            err = 1
        if not self.m_fa:
            sys.stderr.write("Missing argument: --m-fa\n")
            err = 1
        if not self.m_input:
            sys.stderr.write("Missing argument: --m-input\n")
            err = 1
        if not self.c_vcf:
            sys.stderr.write("Missing argument: --c-vcf\n")
            err = 1
        if not self.c_fa:
            sys.stderr.write("Missing argument: --c-fa\n")
            err = 1
        if not self.c_input:
            sys.stderr.write("Missing argument: --c-input\n")
            err = 1
        if not self.c_output:
            sys.stderr.write("Missing argument: --c-output\n")
            err = 1

        if err:
            sys.stderr.write("\n")

        return not err

        

if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "f-vcf=", "f-fa=", "f-input=", "m-vcf=", "m-fa=", "m-input=", "c-vcf=", "c-fa=", "c-input=", "c-output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    args = ArgContainer()
    for o, a in opts:
        if a and a[0] == "=":
            a = a[1:]

     
        if o == "--f-vcf":
            args.f_vcf = a
        elif o == "--f-fa":
            args.f_fa = a
        elif o == "--f-input":
            args.f_input = a
        elif o == "--m-vcf":
            args.m_vcf = a
        elif o == "--m-fa":
            args.m_fa = a
        elif o == "--m-input":
            args.m_input = a
        elif o == "--c-vcf":
            args.c_vcf = a
        elif o == "--c-fa":
            args.c_fa = a
        elif o == "--c-input":
            args.c_input = a
        elif o == "--c-output":
            args.c_output = a
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
