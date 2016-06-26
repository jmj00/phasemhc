import argparse, gzip, sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''
Map variant positions using pairwise alignment of parent-offspring pair and variant info file.
''')

parser.add_argument('-v, --varinfo', type=str, action = 'store', dest = 'info',
                    help='variant info file')
parser.add_argument('-i, --infile', type=str, action = 'store', dest = 'infasta',
                    help='Input fasta file')
parser.add_argument('-o, --outfile', type=str, action = 'store', dest = 'outfasta',
                    help='Output fasta file')
parser.add_argument('-p, --parent', type=str, action = 'store', dest = 'parent',
                    help='parent ("P" for father (paternal), "M" for mother (maternal)')
parser.add_argument('-t, --trio', type=str, action = 'store', dest = 'trio',
                    help='trio')
parser.add_argument('--verbose', action='store_true', help="print verbose output")
args = parser.parse_args()



def fileOpener(f):
    if f.endswith('gz'):
        return gzip.open(f)
    else:
        return(open(f))

def make_record(string, trio, parent):
    """Returns a new SeqRecord with the reverse complement sequence."""
    return SeqRecord(seq = Seq(string), \
                 id = trio + "_" + parent, \
                 description = "")

SEQS = [str(seq_record.seq) for seq_record in SeqIO.parse(fileOpener(args.infasta), 'fasta')]
# save seqs as strings s1 (transmitted) and s2 (non-transmitted) and ensure uppercase
s1 = SEQS[0].upper()
s2 = SEQS[1].upper()
# create list of characters in sequence
s1L = list(s1)
s2L = list(s2)


if args.verbose:
    "Reading alignments..."
    print >>sys.stderr, s1[0:25],"..."
    print >>sys.stderr, s2[0:25],"..."

#variant info file contains header of format:
#chrom pos var len parent method LR P_prob M_prob mtype cqual nreads n_j child_a1 child_a2 father_a1 father_a2 mother_a1 mother_a2

def get_var_info(variant_info_file, parent):
    # Set counters
    var_count = 0
    homo_count = 0
    het_count = 0
    no_mendelian_violation_count = 0
    alt_denovo_count = 0
    ref_denovo_count = 0
    low_parent_GQ_count = 0

    if args.verbose:
        print >>sys.stderr, "Retrieving variant information..."
    var_positions = []
    variants = []

    parent = parent.upper()
    with open(variant_info_file) as f:
        for l in f:
            if not l.startswith('chrom'):
                curr_line = l.strip().split() # read line
                position = curr_line[1] # get variant position in transmitted haplotype
                var = curr_line[2] # get variant in transmitted haplotype
                var_parent = curr_line[4] # get the phasing info ('P' (paternal),'M'(maternal) or 'U' (unassigned))
                mtype = curr_line[9]

                # check whether variant was transmitted from the parent examined in the child-parent pair if parent assigned
                if var_parent == parent.upper():
                    # increment variant count
                    var_count += 1
                    # set genotype according to parent examined
                    if parent == 'P':
                        # check for mendelian violation:
                        if mtype == 'no_mendelian_violation':
                            # increment counter
                            no_mendelian_violation_count += 1
                            # set genotype
                            p_gt = [curr_line[15],curr_line[16]]

                        # set parent genotype so that the first entry is the transmitted allele, the second entry is the non-transmitted allele
                            if var == p_gt[0]:
                                pass
                            elif var == p_gt[1]:
                                p_gt.reverse()
                        else:
                            # Check mendelian error type and increment counter
                            if mtype == 'alt_denovo':
                                alt_denovo_count += 1
                            elif mtype == 'ref_denovo':
                                ref_denovo_count += 1
                            elif mtype == 'low_parent_GQ':
                                low_parent_GQ_count += 1
                            if args.verbose:
                                print >> sys.stderr,"Found mendelian violation. Violation type: ", mtype

                            # Check whether parent is homozygous
                            if p_gt[0] == p_gt[1]:
                                if args.verbose:
                                    print >>sys.stderr,"Parent is homozygous, but mendelian error was found, so setting genotype to ambiguous character (N) of length equal to parental alleles\n----"
                            else:
                                if args.verbose:
                                    print >> sys.stderr,"Setting genotype to ambiguous character (N) of length equal to parental alleles\n----"
                            p_gt == ['N'*len(curr_line[15]), 'N'*len(curr_line[16])]
                            # set parent genotype so that the first entry is 'N'*(length of the transmitted allele), the second entry is'N'*(length of the non-transmitted allele)
                            if len(var) == len(p_gt[0]):
                                pass
                            elif len(var) == len(p_gt[1]):
                                p_gt.reverse()

                    if parent == 'M':
                        # check for mendelian violation:
                        if mtype == 'no_mendelian_violation':
                            # increment counter
                            no_mendelian_violation_count += 1
                            # set genotype
                            p_gt = [curr_line[17],curr_line[18]]

                        # set parent genotype so that the first entry is the transmitted allele, the second entry is the non-transmitted allele
                            if var == p_gt[0]:
                                pass
                            elif var == p_gt[1]:
                                p_gt.reverse()
                        else:
                            # Check mendelian error type and increment counter
                            if mtype == 'alt_denovo':
                                alt_denovo_count += 1
                            elif mtype == 'ref_denovo':
                                ref_denovo_count += 1
                            elif mtype == 'low_parent_GQ':
                                low_parent_GQ_count += 1

                            if args.verbose:
                                print >> sys.stderr,"Found mendelian violation. Violation type: ", mtype

                            # Check whether parent is homozygous
                            if p_gt[0] == p_gt[1]:
                                if args.verbose:
                                    print >>sys.stderr,"Parent is homozygous, but mendelian error was found, so setting genotype to ambiguous character (N) of length equal to parental alleles\n----"
                            else:
                                if args.verbose:
                                    print >> sys.stderr,"Setting genotype to ambiguous character (N) of length equal to parental alleles\n----"
                            p_gt == ['N'*len(curr_line[17]), 'N'*len(curr_line[18])]
                            # set parent genotype so that the first entry is 'N'*(length of the transmitted allele), the second entry is'N'*(length of the non-transmitted allele)
                            if len(var) == len(p_gt[0]):
                                pass
                            elif len(var) == len(p_gt[1]):
                                p_gt.reverse()

                    # Append parent genotype as tuple to list of variants if p_gt is heterozygote
                    if p_gt[0] == p_gt[1]:
                        # increment counter
                        homo_count += 1
                        if args.verbose:
                            print >>sys.stderr,"Parent is homozygous. Skipping variant\n----"
                        continue
                    else:
                        het_count += 1
                        # Append variant position to list of variant positions
                        var_positions.append(int(position))
                        variants.append(tuple(p_gt))

                # if variant is unassigned, set genotype as ambiguous
                if var_parent == 'U':
                    var_count += 1

                    # Append variant position to list of variant positions
                    var_positions.append(int(position))
                    if not mtype == 'no_mendelian_violation':
                        # increment counters
                        if mtype == 'alt_denovo':
                            alt_denovo_count += 1
                        elif mtype == 'ref_denovo':
                            ref_denovo_count += 1
                        elif mtype == 'low_parent_GQ':
                            low_parent_GQ_count += 1
                        if args.verbose:
                            print >> sys.stderr,"Found mendelian violation at unassigned variant. Violation type: ", mtype, "\n----"
                    else:
                        no_mendelian_violation_count += 1

                    if parent == 'P':
                        p_gt = ['N'*len(curr_line[15]), 'N'*len(curr_line[16])]
                        # set parent genotype so that the first entry is 'N'*(length of the transmitted allele), the second entry is'N'*(length of the non-transmitted allele)
                        if len(var) == len(p_gt[0]):
                            pass
                        elif len(var) == len(p_gt[1]):
                            p_gt.reverse()

                    elif parent == 'M':
                        p_gt = ['N'*len(curr_line[17]), 'N'*len(curr_line[18])]
                        if len(var) == len(p_gt[0]):
                            pass
                        elif len(var) == len(p_gt[1]):
                            p_gt.reverse()

                    # Append parent genotype as tuple to list of variants
                    variants.append(tuple(p_gt))

            else:
                continue
    return var_positions, variants, var_count, homo_count, het_count, no_mendelian_violation_count, alt_denovo_count, ref_denovo_count, low_parent_GQ_count

# Get variant information
var_position, variants, var_count, homo_count, het_count, no_mendelian_violation_count, alt_denovo_count, ref_denovo_count, low_parent_GQ_count = get_var_info(args.info, args.parent)

# Variants should be 'variant that was transmitted, variant that was not transmitted (and therefore the variant that should be put in to the sequence if not already there)'
# e.g.:
# variants = [('A','T'),('C','G'), ('TACGG','TG'),('T','TTTTT')]
# var_position = [1,7,12,18]

# Get ungapped positions in transmitted (s1) and non-transmitted (s2) sequences
if args.verbose:
    print >>sys.stderr, "Getting ungapped positions in sequences..."
ungapped_pos_s1 = [i for i in range(len(s1)) if s1[i] != "-"]
ungapped_pos_s2 = [i for i in range(len(s2)) if s2[i] != "-"]


# Check whether variant should be replaced (check whether variant allele is in the sequence already)
if args.verbose:
    print >>sys.stderr, "Replacing variants..."
for i in range(len(var_position)):

    if variants[i][0] == s2[ungapped_pos_s1[var_position[i]-1]:ungapped_pos_s1[var_position[i]-1]+len(variants[i][0])]:
        # Define length of variant allele and allele present in sequence
        ls = len(variants[i][0])
        lr = len(variants[i][1])
        if args.verbose:
            print >>sys.stderr, "Changing variant at position ", var_position[i], "from ", variants[i][0], "to", variants[i][1]

        if lr == ls:
            s2L[ungapped_pos_s1[var_position[i]-1]] = variants[i][1]
        else:
            s2L[ungapped_pos_s1[var_position[i]-1]:ungapped_pos_s1[var_position[i]-1]+len(variants[i][0])] = ["".join([variants[i][1]])]+[""]*abs(ls-1)

    elif variants[i][1] == s2[ungapped_pos_s1[var_position[i]-1]:ungapped_pos_s1[var_position[i]-1]+len(variants[i][1])]:
        if args.verbose:
            print >>sys.stderr, "Variant at position", var_position[i], "is correct. \nSkipping variant..."
        continue
    else:
        # if non of the variant alleles are found in the sequence, skip to next variant and write warning
        if args.verbose:
            print >>sys.stderr, "Warning. Did not find matching variant string at position: ", var_position[i],"\n"
            # check whether variant is ambiguous and if so replace with ambiguity code N of length of variant
        if 'N' in variants[i][0]:
            ls = len(variants[i][0])
            lr = len(variants[i][1])
            if args.verbose:
                print >>sys.stderr,"Changing variant at position ", var_position[i], "to", variants[i][1]
            if lr == ls:
                s2L[ungapped_pos_s1[var_position[i]-1]] = variants[i][1]
            else:
                s2L[ungapped_pos_s1[var_position[i]-1]:ungapped_pos_s1[var_position[i]-1]+len(variants[i][0])] = ["".join([variants[i][1]])]+[""]*abs(ls-1)
        # If not, skip variant
        else:
            if args.verbose:
                print >> sys.stderr,"Variant is non-ambiguous. Skipping variant..."
                # if variant is non-ambiguous decrement het_count
            het_count -= 1
            continue

non_count = var_count - het_count - homo_count

# Create seq record object
record = make_record("".join(s2L).replace("-",""), args.trio, args.parent+"1")

# Write to outfile
with open(args.outfasta, 'w') as out:
    out.write(record.format('fasta'))


print >>sys.stdout, "# Total number of variants:", var_count
print >>sys.stdout, "# Number of phased homozygous variants: ", homo_count
print >>sys.stdout, "# Number of phased heterozygote variants: ", het_count
print >>sys.stdout, "# Number of non-phased variants: ", non_count
print >>sys.stdout, "# ------------------------------"
print >>sys.stdout, "# Number of variants with no mendelian violation:", no_mendelian_violation_count
print >>sys.stdout, "# Number of variants with mendelian error of type alt_denovo: ", alt_denovo_count
print >>sys.stdout, "# Number of variants with mendelian error of type ref_denovo: ", ref_denovo_count
print >>sys.stdout, "# Number of variants with mendelian error of type low_parent_GQ: ", low_parent_GQ_count
print >>sys.stdout, "###############################"
print >>sys.stdout, "# \n\n"
print >>sys.stdout, "# ------------------------------"
print >>sys.stdout, "# Writing in tabular format:"
print >>sys.stdout, "# ------------------------------"
print >>sys.stdout, "\t".join(["var_count","homo_count","het_count","non_count","no_mendelian_violation_count","alt_denovo_count","ref_denovo_count","low_parent_GQ_count","\n"])
print >>sys.stdout, var_count,"\t",homo_count,"\t",het_count,"\t",non_count,"\t",no_mendelian_violation_count,"\t",alt_denovo_count,"\t",ref_denovo_count,"\t",low_parent_GQ_count,"\n"




