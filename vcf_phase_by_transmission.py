import vcf
import sys
import argparse

parser = argparse.ArgumentParser(description='''
Phase a trio using transmission information
''')

parser.add_argument('vcf', type=str,
                    help='VCF-file of called variants.')
parser.add_argument('father', type=str,
                    help='Father of this family.')
parser.add_argument('mother', type=str,
                    help='Mother of this family')
parser.add_argument('child', type=str,
                    help='Child of this family')
parser.add_argument('--verbose', action='store_true', help="print progress output")
parser.add_argument('--region', 
                    default="", action='store', type=str,
                    help='Analyse only the specified region. Either a chromosome: "chr1" or an interval: "chr1:1000-2000"')

args = parser.parse_args()

if args.region == "":
    reader = vcf.Reader(filename=args.vcf)
elif ':' in args.region:
    chrom = args.region.split(":")[0]
    pos1, pos2 = (int(x) for x in args.region.split(":")[1].split('-'))
    reader = vcf.Reader(filename=args.vcf).fetch(chrom, pos1, pos2)   
else:
    reader = vcf.Reader(filename=args.vcf).fetch(args.region)

def check_and_fix(a1, a2, ref, alt):
    if a1 != a2 and a1 != ref:
        assert a2 == ref and a1 == alt
        return ref, alt
    else:
        return a1, a2

for record in reader:
    try:
        if len(record.ALT) != 1:
            print >> sys.stderr, record, "is multiallelic"
            continue
        gt_father = record.genotype(args.father)
        gt_mother = record.genotype(args.mother)
        gt_child = record.genotype(args.child)
        if gt_father.gt_bases is None or gt_mother.gt_bases is None or gt_child.gt_bases is None:
            print >> sys.stderr, record, "one of the genotypes are missing"
            continue            
        father_a1, father_a2 =  gt_father.gt_bases.split('/')
        mother_a1, mother_a2 = gt_mother.gt_bases.split('/')
        child_a1, child_a2 = gt_child.gt_bases.split('/')
        ref = record.REF
        alt = record.ALT[0]
        child_a1, child_a2 = check_and_fix(child_a1, child_a2, ref, alt)
        father_a1, father_a2 = check_and_fix(father_a1, father_a2, ref, alt)
        mother_a1, mother_a2 = check_and_fix(mother_a1, mother_a2, ref, alt)
        print record.CHROM, record.POS, ref, alt,
        if not ((child_a1 in [father_a1, father_a2] and child_a2 in [mother_a1, mother_a2]) or
                (child_a2 in [father_a1, father_a2] and child_a1 in [mother_a1, mother_a2])):
            print 'mendelian_violation',
            print '%s/%s %s/%s %s/%s' %(child_a1, child_a2,
                                        father_a1, father_a2,
                                        mother_a1, mother_a2)
        else:
            print 'no_mendelian_violation',
            if child_a1 == child_a2:
                if child_a1 == ref:
                    print '%s|%s %s|%s %s|%s' %(child_a1, child_a2,
                                                father_a1, father_a2,
                                                mother_a1, mother_a2)
                else:
                    print '%s|%s %s|%s %s|%s' %(child_a1, child_a2,
                                                father_a2, father_a1,
                                                mother_a2, mother_a1)                
            elif (gt_father.gt_type < gt_mother.gt_type):
                print '%s|%s %s|%s %s|%s' %(ref, alt,
                                            father_a1, father_a2,
                                            mother_a2, mother_a1)
            elif (gt_mother.gt_type < gt_father.gt_type):
                print '%s|%s %s|%s %s|%s' %(alt, ref,
                                            father_a2, father_a1,
                                            mother_a1, mother_a2)
            else:
                print '%s/%s %s/%s %s/%s' %(child_a1, child_a2,
                                            father_a1, father_a2,
                                            mother_a1, mother_a2)
    except :
        print >> sys.stderr, record, "Unexpected error,", sys.exc_info()
        continue

