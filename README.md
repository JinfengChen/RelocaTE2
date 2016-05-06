# RelocaTE2: a high resolution mapping tool for transposable elements polymorphisms in large population data

## Introduction

## Installation
+ System requirements
  - Linux/Unix platform
  - Short read aligner: BLAT (v35), bowtie (v2.2.6), bwa (v0.6.2)
  - Python (v2.7.5) and pysam pakcage ()
  - Perl (v5.20.2)
  - seqtk (v1.0)
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



## RelocaTE2 Command Line Options

## RelocaTE2 input

## RelocaTE2 output

## 
