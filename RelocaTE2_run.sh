#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l mem=10gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR

relocate=/rhome/cjinfeng/BigData/00.RD/RelocaTE2/scripts/relocaTE.py
#repeat=/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa
repeat=/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/Rice.TE.short.unique.fa
genome=/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa
ref_te=/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU_r7.fa.RepeatMasker.out
fq_d=/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Landrace/fastq/HEG4
bam=/rhome/cjinfeng/BigData/01.Rice_genomes/HEG4/00.Bam/HEG4_MSU7_BWA/HEG4_2.3.MSU7_BWA.bam
outdir=/rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Real_Data/Landrace/RelocaTEi/HEG4_2
aligner=blat
size=500

start=`date +%s`

python $relocate --te_fasta $repeat --genome_fasta $genome --fq_dir $fq_d --bam $bam --outdir $outdir --reference_ins $ref_te --sample $strain --size $size --step 1234567 --mismatch 2 --run --cpu $PBS_NP --aligner $aligner --split --verbose 3

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
