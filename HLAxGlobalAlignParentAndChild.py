#!/usr/bin/env python
"""

Read fasta input from parent and child.

Read the phased information and filter out any phased entries with positions crossing

Cut the fasta input up, by using the phased information, which provides a position for the child and the parent.

For every section, run global alignment and output the results

Author: Rune Moellegaard Friborg <runef@birc.au.dk>
"""

VCF_OUTPUT_TEMPLATE=""

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

VERSION="0.07"
UPDATED="2015-08-13"
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
        
        self.val = records[0].seq
        self.id = records[0].id
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


class PhasedPositions():
    """
    Read the input of format:
    #CHILDFATHERMOTHER
    3343   6817    3170

    And provide filtering methods
    """
    def __init__(self, filename, parent):

        self.content = []
        
        sys.stdout.write("Parsing '"+filename+"'...")
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
                c0_pos, c1_pos, f0_pos, _ , m0_pos, _ = l
                
                # Typecast

                if parent == 'f':
                    if f0_pos != 'None':
                        p_pos = int(f0_pos)
                        c_pos = int(c0_pos)
                        self.content.append((c_pos, p_pos))
                else:
                    if m0_pos != 'None':
                        p_pos = int(m0_pos)
                        c_pos = int(c1_pos)
                        self.content.append((c_pos, p_pos))

                
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)            

        sys.stdout.write("done\n")
        handle.close()
    
    def clean(self):
        """
        Removes crossed phased positions 
        """
    
        L = []
        
        last = 0
        for c_pos, p_pos in self.content:
            if p_pos > last:
                L.append((c_pos, p_pos))
            else:
                L.pop()
                print(str(last)+ " > " +str(p_pos))

            last = L[-1][1]

        sys.stdout.write("Removed "+str(len(self.content)-len(L))+" cross-phased positions\n")

        self.content = L

class Alignment():
    def __init__(self, fasta, vcf, positions, parent, gapopen, gapextend, max_mem):

        self.vcf = vcf.val
        self.parent = parent
        self.gapopen = gapopen
        self.gapextend = gapextend
        self.max_mem_mb = max_mem*1024
        self.fasta_id = (str(fasta[0].id), str(fasta[1].id))
        self.fasta = (str(fasta[0].val), str(fasta[1].val))
        self.segments = positions


    def merge(self, seq_child, seq_parent):
        
        merged = list(seq_child)
        if (len(seq_child) == len(seq_parent)):
            # merge
            for pos in xrange(len(seq_child)):
                if seq_child[pos] == 'N' or seq_child[pos] == '-':
                    merged[pos] = seq_parent[pos]
        else:
            # Skip beginning and end sequence, which is different in length and have not been aligned.
            pass
        
        return "".join(merged)

    def write_vcf_range(self, vcf_writer, alignment, pos_start, vcf_offset):

        # Write vcf with corrected positions
        pos = pos_start
        for c in alignment:
            if c == '-':
                vcf_offset += 1
            else:
                if self.vcf.has_key(pos):
                    record = self.vcf[pos]
                    record.POS = pos + vcf_offset
                    vcf_writer.write_record(record)
                pos += 1
        
        return pos, vcf_offset


    def process(self, output_prefix):
        
        sys.stdout.write("Main processing:\n")

        vcf_template = vcf.Reader(filename=VCF_OUTPUT_TEMPLATE)

        sys.stdout.write("Running global alignment on "+str(len(self.segments))+" segments\n")
        start = self.segments[0]

        ohandle_fasta_child = open(output_prefix + '.aligned_against_'+self.parent+'.c.fa', "w")
        ohandle_fasta_parent = open(output_prefix + '.aligned_against_c.'+self.parent+'.fa', "w")
        ohandle_merged_parent_into_child = open(output_prefix + '.consensus_of_'+self.parent+'.c.fa', "w")

        ohandle_vcf = open(output_prefix + '.consensus_of_'+self.parent+'.c.vcf', "w")
        vcf_writer = vcf.Writer(ohandle_vcf, vcf_template)

        sys.stdout.write("\tWriting " + output_prefix + '.aligned_against_'+self.parent+".c.fa...\n")
        sys.stdout.write("\tWriting " + output_prefix + '.aligned_against_c.'+self.parent+".fa...\n")
        sys.stdout.write("\tWriting " + output_prefix + '.consensus_of_'+self.parent+'.c.fa...\n')
        sys.stdout.write("\tWriting " + output_prefix + '.consensus_of_'+self.parent+'.c.vcf...\n')

        # Add header
        ohandle_fasta_child.write(">"+self.fasta_id[0]+"\n")
        ohandle_fasta_parent.write(">"+self.fasta_id[1]+"\n")
        ohandle_merged_parent_into_child.write(">"+self.fasta_id[0]+"\n")

        # Add fasta before first segment
        ohandle_fasta_child.write(self.fasta[0][:start[0]])
        ohandle_fasta_parent.write(self.fasta[1][:start[1]])

        ohandle_merged_parent_into_child.write(self.merge(self.fasta[0][:start[0]], self.fasta[1][:start[1]]))

        # Fix VCF offset
        pos_start = 0
        vcf_offset = 0

        pos_start, vcf_offset = self.write_vcf_range(vcf_writer, self.fasta[0][:start[0]], pos_start, vcf_offset)

        for next_phased_pos in self.segments[1:]:
            child = self.fasta[0][start[0]:next_phased_pos[0]]
            parent = self.fasta[1][start[1]:next_phased_pos[1]]
            
            sys.stdout.write("\n"+str(start[0]))
            align_child, align_parent = self.EMBOSSnwalign(child, parent)
            
            i = 0
            j = 0
            for seq_child,seq_parent in zip(align_child, align_parent):
                i += len(seq_child)
                j += len(seq_parent)
                
                ohandle_fasta_child.write(seq_child)
                ohandle_fasta_parent.write(seq_parent)
                
                ohandle_merged_parent_into_child.write(self.merge(seq_child, seq_parent))
                
                pos_start, vcf_offset = self.write_vcf_range(vcf_writer, seq_child, pos_start, vcf_offset)

            sys.stdout.write("Wrote "+str(i)+"|"+str(j)+"\n")
                
            start = next_phased_pos
        
        # Add fasta after last segment
        ohandle_fasta_child.write(self.fasta[0][start[0]:])
        ohandle_fasta_parent.write(self.fasta[1][start[1]:])
        ohandle_merged_parent_into_child.write(self.merge(self.fasta[0][start[0]:],self.fasta[1][start[1]:]))

        pos_start, vcf_offset = self.write_vcf_range(vcf_writer, self.fasta[0][start[0]:], pos_start, vcf_offset)

       
        ohandle_fasta_child.close()
        ohandle_fasta_parent.close()
        ohandle_merged_parent_into_child.close()
        ohandle_vcf.close()


        sys.stdout.write("Main processing completed.\n")        

    def EMBOSSnwalign(self, child, parent, force_stretcher=False):

        estimated_mem = ((14*len(child)*len(parent))/1024)/1024
        sys.stdout.write("..+"+str(len(child))+"|"+str(len(parent)))            

        aSeq = open("/dev/shm/"+PID+"_a", "w")
        bSeq = open("/dev/shm/"+PID+"_b", "w")
        cSeq = open("/dev/shm/"+PID+"_c", "w")
        
        aSeq.write(child)
        aSeq.close()
        bSeq.write(parent)
        bSeq.close()
        cSeq.close()

        if (estimated_mem < self.max_mem_mb and not force_stretcher):
            sys.stdout.write(" mem using needle: " +str(estimated_mem)+"MB\n")
            sys.stdout.flush()

            p = subprocess.Popen(["needle", "-asequence", aSeq.name, "-bsequence", bSeq.name, "-gapopen", str(self.gapopen), "-gapextend", str(self.gapextend), "-datafile", "NUC.4.4", "-outfile",cSeq.name])

        else:
            sys.stdout.write(" using stretcher!! Increase --max-mem to "+str((estimated_mem/1024)+1)+"G to use needle\n")
            sys.stdout.flush()

            p = subprocess.Popen(["stretcher", "-asequence", aSeq.name, "-bsequence", bSeq.name, "-gapopen", str(self.gapopen), "-gapextend", str(self.gapextend), "-datafile", "NUC.4.4", "-outfile",cSeq.name])

        p.communicate(None)

        if p.returncode == 0:
            os.unlink(aSeq.name)
            os.unlink(bSeq.name)

            # Get result
            cSeq = open(cSeq.name, "r")
        
            # Skip first 25 lines
            for i in xrange(25):
                cSeq.next()
            
            sys.stdout.write(cSeq.next())
            sys.stdout.write(cSeq.next())
            sys.stdout.write(cSeq.next())
        
            # Skip 5 lines
            for i in xrange(4):
                cSeq.next()

            result_child = []
            result_parent = []
        
        
            try:
                if (estimated_mem < self.max_mem_mb and not force_stretcher):
                    child_done = False
                    parent_done = False
                    while (not child_done and not parent_done):
                        line = cSeq.next()
                        seq_child = line[21:71]
                        cSeq.next()
                        line = cSeq.next()
                        seq_parent = line[21:71]
                        cSeq.next()
            
                        if ' ' in seq_child:
                            result_child.append(seq_child[:seq_child.index(' ')])
                            child_done = True
                        else:
                            result_child.append(seq_child)
            

                        if ' ' in seq_parent:
                            result_parent.append(seq_parent[:seq_parent.index(' ')])
                            parent_done = True
                        else:
                            result_parent.append(seq_parent)
                else:
                    
                    cSeq.next()

                    child_done = False
                    parent_done = False
                    while (not child_done and not parent_done):
                        line = cSeq.next()
                        seq_child = line[7:57]
                        cSeq.next()
                        line = cSeq.next()
                        seq_parent = line[7:57]
                        cSeq.next()
                        cSeq.next()
                        cSeq.next()

                        if '\n' in seq_child:
                            result_child.append(seq_child[:seq_child.index('\n')])
                            child_done = True
                        else:
                            result_child.append(seq_child)
            

                        if '\n' in seq_parent:
                            result_parent.append(seq_parent[:seq_parent.index('\n')])
                            parent_done = True
                        else:
                            result_parent.append(seq_parent)
                                                
            except StopIteration:
                pass

            os.unlink(cSeq.name)
        else:
            if (estimated_mem < self.max_mem_mb and not force_stretcher):
                sys.stderr.write("Warning! EMBOSS needle failed with returncode: "+ str(p.returncode) + "\n")
                sys.stderr.write("Trying to run with EMBOSS stretcher")
                result_child, result_parent = self.EMBOSSnwalign(child, parent, force_stretcher=True)
            else:
                sys.stderr.write("Error! EMBOSS stretcher failed with returncode: "+ str(p.returncode) + "\n")
                sys.exit(1)

        return result_child, result_parent


##############################################################
####################### Main #################################

def main(args):
    global VCF_OUTPUT_TEMPLATE
    VCF_OUTPUT_TEMPLATE = args.c_vcf

    # Read Fasta files (child and parent)
    faSequence = []
    for filename in [args.c_fa, args.parent_fa]:
        faSequence.append(Fasta(filename))

    # Read VCF file for child
    childVCFRecords = VCF(args.c_vcf)

    phasePos = PhasedPositions(args.phase_pos, args.parent)
    phasePos.clean()

    align = Alignment(faSequence, childVCFRecords, phasePos.content, args.parent, gapopen=args.gapopen, gapextend=args.gapextend, max_mem=args.max_mem)
    align.process(output_prefix=args.c_output_prefix)
    
    
##############################################################
######################### Help ###############################

def usage():
    print("""HLAxGlobalAlignParentAndChild version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxAssembleFastaFromPhasing \
    --parent-vcf=<file> --parent-fa=<fasta file from previous HLAx stage> --parent=m|f \
    --c-vcf=<file> --c-fa=<fasta file from previous HLAx stage> \
    --phase-pos=<output from previous HLAx stage> --c-output-prefix=<file prefix> \
    --gapopen=6 --gapextend=2
    --max-mem=1024 (gigabyte)
""")


class ArgContainer():
    def __init__(self):
        self.parent_vcf    = ""
        self.parent_fa     = ""
        self.parent        = None
        self.c_vcf    = ""
        self.c_fa     = ""
        self.phase_pos  = ""
        self.c_output_prefix  = ""
        self.gapopen = 6
        self.gapextend = 2
        self.max_mem  = 1024

    def ok(self):
        err = 0
        if not self.parent_vcf:
            sys.stderr.write("Missing argument: --parent-vcf\n")
            err = 1
        if not self.parent_fa:
            sys.stderr.write("Missing argument: --parent-fa\n")
            err = 1
        if not self.parent:
            sys.stderr.write("Missing argument: --parent\n")
            err = 1
        if not self.c_vcf:
            sys.stderr.write("Missing argument: --c-vcf\n")
            err = 1
        if not self.c_fa:
            sys.stderr.write("Missing argument: --c-fa\n")
            err = 1
        if not self.phase_pos:
            sys.stderr.write("Missing argument: --c-phase-pos\n")
            err = 1
        if not self.c_output_prefix:
            sys.stderr.write("Missing argument: --c-output-prefix\n")
            err = 1

        if err:
            sys.stderr.write("\n")

        return not err

        

if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "parent-vcf=", "parent-fa=", "parent=", "c-vcf=", "c-fa=", "phase-pos=", "c-output-prefix=", "gapopen=", "gapextend=", "max-mem="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    args = ArgContainer()
    for o, a in opts:
        if a and a[0] == "=":
            a = a[1:]

     
        if o == "--parent-vcf":
            args.parent_vcf = a
        elif o == "--parent-fa":
            args.parent_fa = a
        elif o == "--parent":
            args.parent = a
        elif o == "--c-vcf":
            args.c_vcf = a
        elif o == "--c-fa":
            args.c_fa = a
        elif o == "--phase-pos":
            args.phase_pos = a
        elif o == "--c-output-prefix":
            args.c_output_prefix = a
        elif o == "--gapopen":
            args.gapopen = int(a)
        elif o == "--gapextend":
            args.gapextend = int(a)
        elif o == "--max-mem":
            args.max_mem = int(a)
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
