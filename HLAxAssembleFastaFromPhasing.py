#!/usr/bin/env python
"""

Read input of format:
[root@fe1 HLAassemble]# head 1089.phase.info.txt 
1089 1716 AAAAAAA AAAAAA no_mendelian_violation AAAAAAA/AAAAAA AAAAAAA/AAAAAA AAAAAAA/AAAAAA
1089 3343 C T no_mendelian_violation C|C C|T C|C
1089 9494 A AA no_mendelian_violation A/AA A/AA A/AA
1089 10862 C T no_mendelian_violation C|C C|T C|C
1089 12030 T C no_mendelian_violation T|T T|C T|T
1089 13524 C A no_mendelian_violation C|C C|A C|C

Read the matching VCF and Fasta for each of the individuals in the trio.

Based on the phasing input, output six new fasta files with corresponding VCF files

Author: Rune Moellegaard Friborg <runef@birc.au.dk>
"""

VCF_OUTPUT_TEMPLATE=""

import getopt
import sys
import gzip
from Bio import SeqIO
import vcf
import itertools

##############################################################
####################### Configuration ########################

VERSION="0.05"
UPDATED="2015-07-07"

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

class Assembler():
    """
    Descripe the process.
    """
    def __init__(self, vcf, seq):
        self.vcf = []
        for x in vcf:
            self.vcf.append(x.val)

        self.seq = []
        for x in seq:
            self.seq.append(str(x.val))

        self.phaseinfo = []
        self.phasing = ([], [], [])


    def parseVCFandPhasing(self, gt_vcf_filename, phaseinfo_filename):
        
        filename = phaseinfo_filename
        sys.stdout.write("Parsing '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        # Parse lines
        i = 0
        for line in handle:
            try:
                i += 1

                if line == "\n":
                    continue

                l = line.strip().split(' ')
                _, pos, ref, alt, _, c_gt, f_gt, m_gt = l
                
                # Typecast
                pos = int(pos)

            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)            

            self.phaseinfo.append((pos, ref, alt, c_gt, f_gt, m_gt))
        sys.stdout.write("done\n")
        handle.close()

        filename = gt_vcf_filename
        sys.stdout.write("Parsing '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        # Get phaseinfo generator
        phaseinfo = iter(self.phaseinfo)
        p = next(phaseinfo)
        
        # New structure
        phasedict_c = {}
        phasedict_f = {}
        phasedict_m = {}

        duplicates = {}

        # Parse lines and match position with phaseinfo
        i = 0
        ok = False
        for line in handle:
            try:
                i += 1

                if line == "\n" or line[0] == "#":
                    continue

                l = line.strip().split('\t')
                chrom, pos, _, ref, alt, _, _, _, _, f, m, c = l
                
                # Typecast
                pos = int(pos)
                f_v, f_pos = f.split(':')
                f_pos = int(f_pos)
                m_v, m_pos = m.split(':')
                m_pos = int(m_pos)
                c_v, c_pos = c.split(':')
                c_pos = int(c_pos)
                
            
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)                        
                
            if p[0] == pos:
                new_entry = (p, (chrom, ref, alt), (c_pos, f_pos, m_pos), (c_v, f_v, m_v))
                if phasedict_c.has_key(c_pos) or phasedict_f.has_key(f_pos) or phasedict_m.has_key(m_pos):
                    if phasedict_c.has_key(c_pos):
                        sys.stderr.write("WARNING! Removing multiple phasing to same position:\n")                
                        sys.stderr.write("Child: "+ str(phasedict_c[c_pos]) + "\n")
                        duplicates[phasedict_c[c_pos]] = 1

                    if phasedict_f.has_key(f_pos):
                        sys.stderr.write("WARNING! Removing multiple phasing to same position:\n")                
                        sys.stderr.write("Father: "+ str(phasedict_f[f_pos]) + "\n")
                        duplicates[phasedict_f[f_pos]] = 1

                    if phasedict_m.has_key(m_pos):
                        sys.stderr.write("WARNING! Removing multiple phasing to same position:\n")                
                        sys.stderr.write("Mother: "+ str(phasedict_m[m_pos]) + "\n")
                        duplicates[phasedict_m[m_pos]] = 1
                    
                    duplicates[new_entry] = 1
                else:
                    phasedict_c[c_pos] = new_entry
                    phasedict_f[f_pos] = new_entry
                    phasedict_m[m_pos] = new_entry
                try:
                    p = next(phaseinfo)
                except StopIteration:
                    ok = True
                    break

            elif p[0] < pos:
                sys.stderr.write("Error! VCF information missing for phase info position: "+ str(p[0])+"!\n")
                sys.exit(1)
            else:
                # skip
                continue
            
        if (not ok):
            sys.stderr.write("Error! Missing VCF values for phase info!\n")
            sys.exit(1)

        # Sort phasedict by position and save to self.phasing
        i = 0
        for D in [phasedict_c, phasedict_f, phasedict_m]:            
            keys = D.keys()
            keys.sort()
            for k in keys:
                # check duplicate
                if not duplicates.has_key(D[k]):
                    self.phasing[i].append(D[k])
            i += 1
                                
        sys.stdout.write("done\n")
        handle.close()

    def process(self, output_prefix):

        sys.stdout.write("Main processing:\n")

        vcf_template = vcf.Reader(filename=VCF_OUTPUT_TEMPLATE)

        # Prepare phase_out_pos structure to record the new positions for the six fasta files
        phase_out_pos = {}

        # Open output files
        for index, f in [(0, '.c'),(1, '.f'),(2, '.m')]:
            seq_len = len(self.seq[index])
            for v in [0, 1]:
                ohandle_fasta = open(output_prefix + f + str(v)+ '.fa', "w")
                ohandle_vcf = open(output_prefix + f + str(v) + '.vcf', "w")
                vcf_writer = vcf.Writer(ohandle_vcf, vcf_template)

                sys.stdout.write("\tWriting " + output_prefix + f + str(v) + ".[fa|vcf]...")

                ohandle_fasta.write(">"+self.phasing[index][0][1][0]+f+str(v)+"\n")

                vcf_offset = 0
                region_start = 0
                region_end = 0        
        
                for entry in self.phasing[index]:
                    # Entry format:
                    # (    [1716, 'AAAAAAA', 'AAAAAA', 'AAAAAAA/AAAAAA', 'AAAAAAA/AAAAAA', 'AAAAAAA/AAAAAA'],
                    #      ('1089', 'AAAAAAA', 'AAAAAA'),
                    #      (1716, 5181, 1541),
                    #      ('0/1', '0/1', '0/1')
                    # )                    
                    
                    region_end = entry[2][index]

                    if region_end > seq_len:
                        # variant is after fasta entry, skipping entry!
                        sys.stderr.write("\nWARNING: Ran out of fasta data before finishing all phased entries:\n"+str(entry))
                        sys.stderr.write("\n\n")
                        break

                    # Write prefix to phased variant
                    ohandle_fasta.write(self.seq[index][region_start:region_end])

                    # Get phase value
                    phase_val = entry[0][3+index]

                    # Check for phasing
                    if ('|' in phase_val):

                        # Find reference length to replace
                        ref_index, _ = entry[3][index].split('/')
                        ref = entry[1][1+int(ref_index)]
                        replace_len = len(ref)
                        
                        # Write variant
                        val = phase_val.split('|')[v]
                        ohandle_fasta.write(val)

                        new_start = region_end + replace_len
                        vcf_offset = vcf_offset + len(val) - replace_len
                        
                        # Add content to phase_out
                        original_c_pos = entry[0][0]
                        if not phase_out_pos.has_key(original_c_pos):
                            phase_out_pos[original_c_pos] = [None, None, None, None, None, None]
                        phase_out_pos[original_c_pos][index*2+v] = region_end+vcf_offset

                    else:
                        # No phasing!
                        # Keep fasta content and no change to vcf offset
                        new_start = region_end                    
                        

                    # Write vcf with corrected positions
                    for pos in xrange(region_start, new_start):
                        if self.vcf[index].has_key(pos):
                            record = self.vcf[index][pos]
                            record.POS = record.POS + vcf_offset
                            vcf_writer.write_record(record)
                            record.POS = record.POS - vcf_offset

                    # Update region start
                    region_start = new_start

                # Add suffix to last phased variant
                s = self.seq[index][region_start:]
                ohandle_fasta.write(s)


                # Update vcf for last entries
                for pos in xrange(region_start, region_start+len(s)):
                        if self.vcf[index].has_key(pos):
                            record = self.vcf[index][pos]
                            record.POS = record.POS + vcf_offset
                            vcf_writer.write_record(record)
                            record.POS = record.POS - vcf_offset

                ohandle_vcf.close()
                ohandle_fasta.write("\n")
                ohandle_fasta.close()
        
                sys.stdout.write("done\n")
        
        # Write phase out list
        phase_out_handle = open(output_prefix + ".phase.pos", "w")
        phase_out_handle.write("#CHILD0\tCHILD1\tFATHER0\tFATHER1\tMOTHER0\tMOTHER1\n")
        
        keys = phase_out_pos.keys()
        keys.sort()

        for original_phase_pos in keys:
            phase_out_handle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(*phase_out_pos[original_phase_pos]))

        #for c0,c1,f0,f1,m0,m1 in itertools.izip_longest(*phase_out_list):
        #    phase_out_handle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(c0,c1,f0,f1,m0,m1))
        
        phase_out_handle.close()
        sys.stdout.write("Main processing completed.\n")        
        
##############################################################
####################### Main #################################

def main(args):
    global VCF_OUTPUT_TEMPLATE
    VCF_OUTPUT_TEMPLATE = args.c_vcf

    # Read vcf file (child, father, mother)
    vcfRecords = []
    for filename in [args.c_vcf, args.f_vcf, args.m_vcf]:
        vcfRecords.append(VCF(filename))

    # Read Fasta files (child, father, mother)
    faSequence = []
    for filename in [args.c_fa, args.f_fa, args.m_fa]:
        faSequence.append(Fasta(filename))

    # The main object for storing and assembling the haplotypes.
    assembler = Assembler(vcfRecords, faSequence)

    # Parse GT VCF input and phasing
    assembler.parseVCFandPhasing(args.c_gt_vcf, args.c_input)

    # Write haplotype fasta files and corrected VCF files
    assembler.process(output_prefix=args.c_output_prefix)


##############################################################
######################### Help ###############################

def usage():
    print("""HLAxAssembleFastaFromPhasing version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxAssembleFastaFromPhasing --f-vcf=<file> --f-fa=<fasta file> \
    --m-vcf=<file> --m-fa=<fasta file> \
    --c-vcf=<file> --c-fa=<fasta file> --c-gt-vcf=<output from previous HLAx stage> --c-input=<input file> --c-output-prefix=<file prefix>
""")


class ArgContainer():
    def __init__(self):
        self.f_vcf    = ""
        self.f_fa     = ""
        self.m_vcf    = ""
        self.m_fa     = ""
        self.c_vcf    = ""
        self.c_fa     = ""
        self.c_gt_vcf  = ""
        self.c_input  = ""
        self.c_output_prefix  = ""

    def ok(self):
        err = 0
        if not self.f_vcf:
            sys.stderr.write("Missing argument: --f-vcf\n")
            err = 1
        if not self.f_fa:
            sys.stderr.write("Missing argument: --f-fa\n")
            err = 1
        if not self.m_vcf:
            sys.stderr.write("Missing argument: --m-vcf\n")
            err = 1
        if not self.m_fa:
            sys.stderr.write("Missing argument: --m-fa\n")
            err = 1
        if not self.c_vcf:
            sys.stderr.write("Missing argument: --c-vcf\n")
            err = 1
        if not self.c_fa:
            sys.stderr.write("Missing argument: --c-fa\n")
            err = 1
        if not self.c_gt_vcf:
            sys.stderr.write("Missing argument: --c-gt-vcf\n")
            err = 1
        if not self.c_input:
            sys.stderr.write("Missing argument: --c-input\n")
            err = 1
        if not self.c_output_prefix:
            sys.stderr.write("Missing argument: --c-output-prefix\n")
            err = 1

        if err:
            sys.stderr.write("\n")

        return not err

        

if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "f-vcf=", "f-fa=", "m-vcf=", "m-fa=", "c-vcf=", "c-fa=", "c-gt-vcf=", "c-input=", "c-output-prefix="])
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
        elif o == "--m-vcf":
            args.m_vcf = a
        elif o == "--m-fa":
            args.m_fa = a
        elif o == "--c-vcf":
            args.c_vcf = a
        elif o == "--c-fa":
            args.c_fa = a
        elif o == "--c-gt-vcf":
            args.c_gt_vcf = a
        elif o == "--c-input":
            args.c_input = a
        elif o == "--c-output-prefix":
            args.c_output_prefix = a
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
