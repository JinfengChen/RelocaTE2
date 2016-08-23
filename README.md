# RelocaTE2: a high resolution mapping tool for transposable elements polymorphisms in large population data

## Introduction
RelocaTE2 is an improved version of RelocaTE ([Robb et al., 2013](http://www.g3journal.org/content/3/6/949.long)). RelocaTE2 is highly sensitive and accurate in mapping transposable elements (TE) polymorphisms at single base pair resolution. RelocaTE2 uses the reads associated with TEs as seeds to cluster the read pairs on chromosomes. It automatically detects the target site duplication (TSD) of a TE insertion from alignments in each cluster, which enable high resolution mapping of TE polymorphisms. Unlike parallel searching of multi-TE elements in RelocaTE, RelocaTE2 searches all TEs in one cycle, which enable us find polymorphisms of thousands of TEs in an individual genome or large populations in reasonable timeframe without losing sensitivity and specificity.

## Installation
+ System requirements
  - Linux/Unix platform
  - Short read aligner: BLAT (v35+), bowtie2 (v2.2.6+), bwa (v0.6.2)
  - Python (v2.7.5+) and pysam package (v0.8.5+)
  - Perl (v5.20.2+)
  - seqtk (v1.0+)

+ Install
```shell
git clone https://github.com/JinfengChen/RelocaTE2.git
cd RelocaTE2
bash install.sh
uz test_data.tar.gz
cd test_data
bash run_test.sh > run_test.sh.log 2>&1 &
```

+ Troubleshooting
  - installation of RelocaTE2 using install.sh will install all the tools and packages required to run RelocaTE2. The script install and link the executables of all tools to RelocaTE2/bin directory and record their paths in RelocaTE2/CONFIG. The main script of RelocaTE2, relocaTE.py, search executables in $PATH, however, overwrite the values in $PATH using executables from RelocaTE2/CONFIG. Users can modify RelocaTE2/CONFIG with tools installed in their system to avoid any problem. 
  - Python module "pysam" is installed to RelocaTE2/pythonlib. By setting PYTHONPATH=RelocaTE2/pythonlib/lib64/python2.7/site-packages, any old versions installed in the system is overwritten temporary and use the new version for RelocaTE2.

## Quick Start Quide
  - index reference genome
```shell
#reference genome
ref=test_data/FLY603.Chr2L.fa
bwa index $ref
```
  - index repeat sequence if using bowtie2 as search engine (default is to use BLAT as search engine)
```shell
#repeat elements
repeat=test_data/pogo.fa
bowtie2-build $repeat $repeat
```
  - run RelocaTE2 to find transposable element insertions
```shell
#repeatmasker results of TE annotation on reference genome
ref_te=test_data/FLY603.Chr2L.fa.RepeatMasker.out
#directory where the input fastq format reads are located
fq_dir=test_data/FLY603.Chr2L.pogo.rep1_reads/
#output directory where RelocaTE2 write temperary and final output
outdir=test_data/FLY603.Chr2L.pogo.rep1_RelocaTE2_outdir
python relocaTE.py --te_fasta $repeat --genome_fasta $ref --fq_dir $fq_dir --outdir $outdir --reference_ins $ref_te
```

## RelocaTE2 Command Line Options

```shell
python relocaTE2.py
usage: relocaTE.py [-h] [-b BAM] [-t TE_FASTA] [-d FQ_DIR] [-g GENOME_FASTA]
                   [-r REFERENCE_INS] [-o OUTDIR] [-s SIZE] [-c CPU]
                   [-1 MATE_1_ID] [-2 MATE_2_ID] [-u UNPAIRED_ID]
                   [--sample SAMPLE] [--aligner ALIGNER]
                   [--len_cut_match LEN_CUT_MATCH]
                   [--len_cut_trim LEN_CUT_TRIM] [--mismatch MISMATCH]
                   [--mismatch_junction MISMATCH_JUNCTION] [--step STEP]
                   [--run] [--split] [-v VERBOSE]

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     Name of BAM file of reads mapped reference genome
  -t TE_FASTA, --te_fasta TE_FASTA
                        Name of fasta sequence of repeat element
  -d FQ_DIR, --fq_dir FQ_DIR
                        Name of directory of input fastq sequence data
  -g GENOME_FASTA, --genome_fasta GENOME_FASTA
                        Name of fasta file of reference genome sequence
  -r REFERENCE_INS, --reference_ins REFERENCE_INS
                        Name of RepeatMasker TE annotation of reference genome
  -o OUTDIR, --outdir OUTDIR
                        Name of output directory where to put temperary and
                        final results
  -s SIZE, --size SIZE  Insert size of sequence library, default = 500
  -c CPU, --cpu CPU     Number of CPUs to use for multiplex, default = 1
  -1 MATE_1_ID, --mate_1_id MATE_1_ID
                        string define paired-end read1, default = "_1"
  -2 MATE_2_ID, --mate_2_id MATE_2_ID
                        string define paired-end read2, default = "_2"
  -u UNPAIRED_ID, --unpaired_id UNPAIRED_ID
                        string defining single-end reads, default = ".unPaired"
  --sample SAMPLE       string defining sample name which will present in output
                        GFF, default = "not_given"
  --aligner ALIGNER     aligner used to map reads to repeat elements,
                        default=blat
  --len_cut_match LEN_CUT_MATCH
                        length cutoff threshold for match between reads and
                        repeat elements. Large value will lead to less
                        sensitive but more accuracy, default = 10
  --len_cut_trim LEN_CUT_TRIM
                        length cutoff threshold for trimed reads after
                        trimming repeat sequence from reads. Large value will
                        lead to less sensitive but more accuracy, default = 10
  --mismatch MISMATCH   Number of mismatches allowed for matches between reads
                        and repeat elements, default = 2
  --mismatch_junction MISMATCH_JUNCTION
                        Number of mismatches allowed for matches between
                        junction reads and repeat elements, default = 2
  --step STEP           Number to control steps of pipeline, default =
                        "1234567"
  --run                 run while this script excute
  --split               split fastq into 1 million reads chunks to run blat/bwa jobs
  -v VERBOSE, --verbose VERBOSE
                        verbose grade to print out information in all scripts:
                        range from 0 to 4, default = 2
```

## RelocaTE2 input
+ Reference genome sequence (ref): multiple sequences of reference genome in fasta format
```shell
cat test_data/FLY603.Chr2L.fa | head -n 4
>Chr2L type=golden_path_region; loc=2L:1..23513712; ID=2L; dbxref=GB:AE014134,GB:AE014134,REFSEQ:NT_033779; MD5=b6a98b7c676bdaa11ec9521ed15aff2b; length=23513712; release=r6.03; species=Dmel;
CGACAATGCACGACAGAGGAAGCAGAACAGATATTTAGATTGCCTCTCATTTTCTCTCCCATATTATAGGGAGAAATATG
ATCGCGTATGCGAGAGTAGTGCCAACATATTGTGCTCTTTGATTTTTTGGCAACCCAAAATGGTGGCGGATGAACGAGAT
GATAATATATTCAAGTTGCCGCTAATCAGAAATAAATTCATTGCAACGTTAAATACAGCACAATATATGATCGCGTATGC
```
+ Repeat sequence (repeat): Consensus sequences of repeat families in fasta format
```shell
cat test_data/pogo.fa | head -n 4
>pogo
CAGTATAATTCGCTTAGCTGCATCGATAGTTAGCTGCATCGGCAAGATAT
CTGCATTATTTTTCCATTTTTTTGTGTGAATAGAAAATTTGTACGAAAAT
TCATACGTTTGCTGCATCGCAGATAACAGCCTTTTTAACTTAAGTGCATC
```
+ Resequencing data (fq_dir): Illumina reads from one strain in fastq format (fastq or fastq.gz). Sequences need to be put in one directory. Paired-end reads need to end with \_1.fastq and \_2.fastq.
```shell
ls test_data/FLY603.Chr2L.pogo.rep1_reads/*.fq
test_data/FLY603.Chr2L.pogo.rep1_reads/FLY603.Chr2L.pogo.rep1_reads_2X_100_500_1.fq
test_data/FLY603.Chr2L.pogo.rep1_reads/FLY603.Chr2L.pogo.rep1_reads_2X_100_500_2.fq

cat test_data/FLY603.Chr2L.pogo.rep1_reads/FLY603.Chr2L.pogo.rep1_reads_2X_100_500_1.fq | head -n 4
@read_500_1/1
TAATTTTCAATCTTACAGAACTGTTGCGTGGCAAATCCTTCAGCTGGAACGATCATGCTCAAGAAGCTTTTGACAATATCAAAGACAAGTTATGCTCTGC
+
GGGGGGFCCFGHEEE:?<8;C:;@<6=@BFGHEEHHHGGEBEEEEFFH?F;<<@HHFECFFHHHGFFFHHHFBE<BCDCDFHFHHHHHGBHHEED3;DHH

cat test_data/FLY603.Chr2L.pogo.rep1_reads/FLY603.Chr2L.pogo.rep1_reads_2X_100_500_2.fq | head -n 4
@read_500_1/2
GTCAACATCCTCGAACGATCGAGACAAAGCGTCCGCAACCACATTTTCTGTTCCATTGCGATGCTCGATTTCGAAGGTATAACCCTGAAGTTTGATAGCC
+
FFHHHHHHHHHGGGHHHHHHHEHFFEHFHHGHHHDHHHGHFFGEHHFHHHFGDC@=@GFGB4DEHHHECE>>8/@ABHHH?DGHBEFEDAEDHFBDBCC>
```
+ RepeatMasker results of TE annotation on reference genome (ref_te): default TE annotation of reference genome used to call TE insertions in reference genome and remove possible false positive non-reference TE insertions.
```shell
cat test_data/FLY603.Chr2L.fa.RepeatMasker.out | head -n 4
43675 0.48 0.02 0.02                Chr2L     47514     52519 (23461193) +               jockey         Unknown       2    5007    (13)    
   239 19.30 0.00 0.00                Chr2L    215151    215207 (23298505) C                  roo         Unknown  (7976)    1116    1060    
   227 26.05 1.56 7.03                Chr2L    239481    239608 (23274104) C                  roo         Unknown  (7911)    1181    1061    
 67696 0.16 0.01 0.46                Chr2L    347941    355383 (23158329) C                blood         Unknown     (0)    7410       1    
```
+ Bam file of reads mapping results on reference genome (bam): characterize heterozygous and homozygous TE insertions when BAM file is provided.

## RelocaTE2 output
+ Structure of output directory

+ TE insertions
  - TE insertions shared between resequenced strain and reference genome: test\_data/FLY603.Chr2L.pogo.rep1\_RelocaTE2\_outdir/repeat/results/ALL.all\_ref\_insert.gff
  - TE insertions only present in resequenced strain: test\_data/FLY603.Chr2L.pogo.rep1\_RelocaTE2\_outdir/repeat/results/ALL.all_nonref_insert.gff
  - TE insertions characterized as heterozygous and homozygous as described in [Robb et al., 2013](http://www.g3journal.org/content/3/6/949.long): test\_data/FLY603.Chr2L.pogo.rep1\_RelocaTE2\_outdir/repeat/results/ALL.all_nonref_insert.characTErized.gff.
+ GFF format used in RelocaTE2

```shell
Chr2L   RelocaTE2       FLY603  65072   65076   .       -       .       ID=repeat_Chr2L_65072_65076;Name=pogo;TSD=AGAAC;Note=Non-reference, not found in reference;Right_junction_reads=3;Left_junction_reads=1;Right_support_reads=4;Left_support_reads=2;
Chr2L   RelocaTE2       FLY603  198322  198326  .       -       .       ID=repeat_Chr2L_198322_198326;Name=pogo;TSD=ATCCA;Note=Non-reference, not found in reference;Right_junction_reads=1;Left_junction_reads=2;Right_support_reads=2;Left_support_reads=5;
Chr2L   RelocaTE2       FLY603  246039  246043  .       -       .       ID=repeat_Chr2L_246039_246043;Name=pogo;TSD=AAAGG;Note=Non-reference, not found in reference;Right_junction_reads=2;Left_junction_reads=3;Right_support_reads=2;Left_support_reads=5;
Chr2L   RelocaTE2       FLY603  423544  423548  .       +       .       ID=repeat_Chr2L_423544_423548;Name=pogo;TSD=GTGCA;Note=Non-reference, not found in reference;Right_junction_reads=1;Left_junction_reads=1;Right_support_reads=3;Left_support_reads=2;
```
Attributes in field 8 of GFF:

ID: unique id of TE insertions, repeat\_"chromosome"\_"start"\_"end"

Name: TE family name of this insertion

TSD: target site duplicate predicted from read alignments

Note: definition of TE insertions, including Non-reference, Reference-Only and Shared.

Right\_junction\_reads: number of reads covering the junction of TE insertion on right side/downstream.

Left\_junction\_reads: number of reads covering the junction of TE insertion on left side/upstream.

Right\_support\_reads: number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on right side/downstream.

Left\_support\_reads: number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on left side/downstream.

## Publications
1. Robb S.M., Lu L., Valencia E., Burnette J.M. 3rd., Okumoto Y., Wessler S.R., Stajich J.E. The use of RelocaTE and unassembled short reads to produce high-resolution snapshots of transposable element generated diversity in rice. G3 2013;3:949-957.
