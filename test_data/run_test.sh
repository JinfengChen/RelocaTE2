#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

current_dir=`pwd`
echo "$current_dir"
relocate=$current_dir/../scripts/relocaTE2.py
repeat=$current_dir/RiceTE.fa
genome=$current_dir/MSU7.Chr3_2M.fa
ref_te=$current_dir/MSU7.Chr3_2M.fa.RepeatMasker.out
fq_d=$current_dir/MSU7.Chr3_2M.ALL_reads
outdir=$current_dir/MSU7.Chr3_2M.ALL_reads_RelocaTE2_outdir
aligner=blat
size=500
strain=rice

start=`date +%s`
export PYTHONPATH=../lib/python2.7/site-packages
python $relocate --te_fasta $repeat --genome_fasta $genome --fq_dir $fq_d --outdir $outdir --reference_ins $ref_te --sample $strain --size $size --step 1234567 --mismatch 2 --cpu 1 --aligner $aligner --verbose 4
#bedtools window -w 10 -a MSU7.Chr3_2M.ALL.gff -b MSU7.Chr3_2M.ALL_reads_RelocaTE2_outdir/repeat/results/ALL.all_nonref_insert.gff | wc -l

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
