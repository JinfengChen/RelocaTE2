# RelocaTE2: a high resolution mapping tool for transposable elements polymorphisms in large population data

## Introduction
RelocaTE2 is an improved version of RelocaTE ([Robb et al., 2013](http://www.g3journal.org/content/3/6/949.long)). RelocaTE2 is highly sensitivity and accuracy for mapping transposoable elements (TE) polymorphisms at single base pair resolution. RelocaTE2 use the reads that associated with TEs as seeds to cluster the read pairs on chromosomes. It automatically detects target site duplication (TSD) of TE insertion from alignments in each cluster, which enable high resolution mapping of TE polymorphisms. Unlike parallel searching of multi-TE elements in RelocaTE, RelocaTE2 search all TEs in one cycle, which enable us find polymorphisms of thousands of TEs in individual genome or large populations in reasonable timeframe without losing sensitivity and specificity.

## Installation
+ System requirements
  - Linux/Unix platform
  - Short read aligner: BLAT (v35+), bowtie (v2.2.6+), bwa (v0.6.2)
  - Python (v2.7.5+) and pysam pakcage (v0.8.5+)
  - Perl (v5.20.2+)
  - seqtk (v1.0+)

+ Install
```shell
git clone https://github.com/JinfengChen/RelocaTE2.git
cd RelocaTE2
#edit fullpath of tools in CONFIG and run testing script
vi CONFIG
uz test_data.tar.gz
cd test_data
bash RelocaTE2_run.sh > RelocaTE2_run.log 2>&1 &
```
+ Troublethrough
  - install and setup pysam in local directory
```shell
git clone https://github.com/pysam-developers/pysam.git
cd pysam
python setup.py install --prefix ~/software/tools/pythonlib
export PYTHONPATH=$PYTHONPATH:~/BigData/software/pythonlib/lib/python2.7/site-packages
```
## Quick Start Quide
  - index referene genome
```shell
#reference genome
ref=test_data/FLY603.Chr2L.fa
bwa index $ref
```
  - index repeat sequence if use bowtie2 as search engine (default is to use BLAT as search engine)
```shell
#repeat elements
repeat=test_data/pogo.fa
bowtie2-build $repeat $repeat
```
  - run RelocaTE2 to find transposable element insertions.
```shell
#repeatmasker results of TE annotation on reference genome
ref_te=test_data/FLY603.Chr2L.fa.RepeatMasker.out
#directory where to put input fastq format reads
fq_dir=test_data/FLY603.Chr2L.pogo.rep1_reads/
#output directory where RelocaTE2 write temperary and final output
outdir=test_data/FLY603.Chr2L.pogo.rep1_RelocaTE2_outdir
python relocaTE.py --te_fasta $repeat --genome_fasta $ref --fq_dir $fq_dir --outdir $outdir --reference_ins $ref_te 
```

## RelocaTE2 Command Line Options

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
+ Resequencing data (fq_dir): Illumina reads of one strain in fastq format (fastq or fastq.gz). Sequence need to be put in one directory. Paired-end reads need to end with \_1.fastq and \_2.fastq.
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
+ Bam file of read mapping results on reference genome (bam): characterize heterozygous and homozygous TE insertions when BAM file is provided.
## RelocaTE2 output

## 
