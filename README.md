# phasemhc
Scripts used in pipeline to phase mhc haplotypes using transmission- and read-backed phasing and exact alignment of variant graph bubbles within the trio. All scripts are provided as is. Below follows a brief description of the individual scripts.

### Align variants or variants bubbles within trio

    usage: alignVar.py [-h] [--father_fasta FATHER_FASTA]
                   [--mother_fasta MOTHER_FASTA] [--child_fasta CHILD_FASTA]
                   [--father_vcf FATHER_VCF] [--mother_vcf MOTHER_VCF]
                   [--child_vcf CHILD_VCF] [--father_out FATHER_OUT]
                   [--mother_out MOTHER_OUT] [--child_out CHILD_OUT]
                   [--verbose]

    Find positions of variants in a trio through exact matching of upstream
    sequences

    optional arguments:
    -h, --help            show this help message and exit
    --father_fasta FATHER_FASTA 
                        fasta file from father
    --mother_fasta MOTHER_FASTA
                        fasta file from mother
    --child_fasta CHILD_FASTA
                         fasta file from child
    --father_vcf FATHER_VCF
                        vcf file from father
    --mother_vcf MOTHER_VCF
                         vcf file from mother
    --child_vcf CHILD_VCF
                         vcf file from child
    --father_out FATHER_OUT
                         output file from father
    --mother_out MOTHER_OUT
                         output file from mother
    --child_out CHILD_OUT
                         output file from child
    --verbose             print verbose output

## Parse multiple alignment files (MAF) files, concatenate and order scaffolds and update variant files accordingly.

    usage: parseMAF.py [-h] [-m MAF] [-s SCAFFOLD_LIST] [-i IPREFIX] [-o OPREFIX]
                       [-c CPREFIX] [--HLA-F HLAF_FASTA] [--KIFC1 KIFC1_FASTA]
                       [--blast-dir BLAST_DIR]
    
    Concatenate scaffold fastas in correct order, by using alignmentblocks from
    MAF file. Adjust VCFs accordingly
    
    optional arguments:
      -h, --help            show this help message and exit
      -m MAF, --maf MAF     maf file of alignment of scaffolds to reference
      -s SCAFFOLD_LIST, --scaffold-list SCAFFOLD_LIST
                            list of scaffolds
      -i IPREFIX, --in-prefix IPREFIX
                            input prefix of fasta and vcf file
      -o OPREFIX, --out-prefix OPREFIX
                            output prefix of fasta and vcf file
      -c CPREFIX, --complete-prefix CPREFIX
                            prefix of comcplete (concatenated) fasta and vcf file
      --HLA-F HLAF_FASTA    HLA-F fasta file
      --KIFC1 KIFC1_FASTA   KIFC1 fasta file
      --blast-dir BLAST_DIR
                            Directory of BLAST binaries
                            
## Create global alignment between variant anchors of parent-offspring sequences to create consensus sequence.

    HLAxGlobalAlignParentAndChild version 0.07 by Rune M. Friborg (updated 2015-08-13)
    Usage:
      HLAxAssembleFastaFromPhasing     --parent-vcf=<file> --parent-fa=<fasta file from previous HLAx stage> --parent=m|f     --c-vcf=<file> --c-fa=<fasta file from previous HLAx stage>     --phase-pos=<output from previous HLAx stage> --c-output-prefix=<file prefix>     --gapopen=6 --gapextend=2
        --max-mem=1024 (gigabyte)
        
## Assemble fasta sequence from phasing info

    HLAxAssembleFastaFromPhasing version 0.05 by Rune M. Friborg (updated 2015-07-07)
    Usage:
      HLAxAssembleFastaFromPhasing --f-vcf=<file> --f-fa=<fasta file>     --m-vcf=<file> --m-fa=<fasta file>     --c-vcf=<file> --c-fa=<fasta file> --c-gt-vcf=<output from previous HLAx stage> --c-input=<input file> --c-output-prefix=<file prefix>
      
## Correct consensus fasta file and vcf-file to update positions according to insertions/deletions of variatnts

    Usage:
      HLAxCorrectConsensusFasta     --fa-input=<fasta file> --var-pos=<file:column> --var-ref-len=<file:column>
        --var-new-seq=<file:column> --var-new-alt=<file:column-offset> --var-vcf-info=<file:column>
        --fa-output=<fasta file> --vcf-output=<vcf file>
        
## Prepare variant info for phasing and rescue variants within trio

        Usage:
          HLAxCorrectConsensusFasta     --fa-input=<fasta file> --var-pos=<file:column> --var-ref-len=<file:column>
            --var-new-seq=<file:column> --var-new-alt=<file:column-offset> --var-vcf-info=<file:column>
            --fa-output=<fasta file> --vcf-output=<vcf file>

## Map variant positions in parent-offspring pair (to map variants in transmitted to non-transmitted haplotypes)

    usage: map_variants.py [-h] [-v, --varinfo INFO] [-i, --infile INFASTA]
                           [-o, --outfile OUTFASTA] [-p, --parent PARENT]
                           [-t, --trio TRIO] [--verbose]
    
    Map variant positions using pairwise alignment of parent-offspring pair and
    variant info file.
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --varinfo INFO    variant info file
      -i, --infile INFASTA  Input fasta file
      -o, --outfile OUTFASTA
                            Output fasta file
      -p, --parent PARENT   parent ("P" for father (paternal), "M" for mother
                            (maternal)
      -t, --trio TRIO       trio
      --verbose             print verbose output

## Phase variants by transmission

    usage: vcf_phase_by_transmission.py [-h] [--verbose] [--region REGION]
                                        vcf father mother child
    
    Phase a trio using transmission information
    
    positional arguments:
      vcf              VCF-file of called variants.
      father           Father of this family.
      mother           Mother of this family
      child            Child of this family
    
    optional arguments:
      -h, --help       show this help message and exit
      --verbose        print progress output
      --region REGION  Analyse only the specified region. Either a chromosome:
                       "chr1" or an interval: "chr1:1000-2000"
                       
## Determine parent of origin (transmission and read-backed)

    usage: vcf_parent_of_origin.py [-h] [--verbose]
                                   [--min-parents-GQ MIN_PARENTS_GQ]
                                   [--min-child-GQ MIN_CHILD_GQ] [--region REGION]
                                   [--max-marker-distance MAX_MARKER_DISTANCE]
                                   [--cpp CPP]
                                   vcf father mother child child_bam
    
    Assign parent of origin to heterozygous variants in a child using read-backed
    phasing
    
    positional arguments:
      vcf                   VCF-file of called variants.
      father                Father of this family.
      mother                Mother of this family
      child                 Child of this family
      child_bam             Bam file with the reads from the child
    
    optional arguments:
      -h, --help            show this help message and exit
      --verbose             print progress output
      --min-parents-GQ MIN_PARENTS_GQ
                            Minimum GQ for variants in parents when assigning
                            obvious phase.
      --min-child-GQ MIN_CHILD_GQ
                            Minimum GQ for variants in child to be concidered.
      --region REGION       Analyse only the specified region. Either a
                            chromosome: "chr1" or an interval: "chr1:1000-2000"
      --max-marker-distance MAX_MARKER_DISTANCE
                            Maximal distance between two het-markers that will be
                            considered.
      --cpp CPP             Specify the child-parent pair consensus sequence
                            examined (cf/cm)
