#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob
import multiprocessing as mp

def usage():
    test="name"
    message='''
RelocaTEi: improved version of RelocaTE for calling transposable element insertions

bam mode:
python relocaTE.py --bam MSU7.Chr4.ALL.rep1_reads_2X_100_500.bam --genome_fasta MSU7.Chr4.fa --te_fasta mping.fa --reference_ins MSU_r7.fa.RepeatMasker.Chr4.out --outdir RelocaTE_output_multiTE_bam

fastq mode:
python relocaTE.py --fq_dirMSU7.Chr4.ALL.rep1_reads_5X_100_500 --genome_fasta MSU7.Chr4.fa --te_fasta mping.fa --reference_ins MSU7.Chr4.fa.RepeatMasker.out --outdir RelocaTE_output_mPing_gz


--cpu: cpu numbers to run multiprocess jobs, default=1
--run: run all steps while excute this script (1) or only generate scripts for all steps (0), default=0

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def writefile(outfile, lines):
    ofile = open(outfile, 'w')
    print >> ofile, lines
    ofile.close()


def shell_runner(cmdline):
    try:
        os.system(cmdline)
    except:
        return 0
    return 1

##run multi process job using pool with limited number of cpu
##cmds is list of shell command, cpu is number of cpu to use
def mp_pool(cmds, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(shell_runner, cmds)
    count= 0
    for x in imap_it:
        print 'job: %s' %(cmds[count])
        print 'status: %s' %(x)
        count += 1

##run job by sequence
def single_run(cmds):
    for cmd in cmds:
        status = shell_runner(cmd)
        print 'job: %s' %(cmd)
        print 'status: %s' %(status)

def existingTE_RM_ALL(top_dir, infile):
    ofile_RM = open('%s/existingTE.bed' %(top_dir), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\s+',line)
                if not unit[0] == '':
                    unit.insert(0, '')
                #print line
                #print unit[5], unit[9], unit[12], unit[13], unit[14]
                if unit[9] == '+':
                    #for i in range(int(unit[6])-2, int(unit[6])+3):
                    #    existingTE_inf[unit[5]]['start'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7], '1', '+')
                    #print unit[10], 'start', unit[6]
                    #for i in range(int(unit[7])-2, int(unit[7])+3):
                    #    existingTE_inf[unit[5]]['end'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '+')
                    #print unit[10], 'end', unit[7]

                    ##if this repeat is a intact element
                    intact = 0
                    if int(unit[12]) == 1 and len(unit[14]) == 3:
                        unit[14] =re.sub(r'\(|\)', '', unit[14])
                        if int(unit[14]) == 0:
                            intact = 1
                    print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])), str(int(unit[7])), unit[11],unit[6],unit[7], intact, '+')
                elif unit[9] == 'C':
                    #for i in range(int(unit[6])-2, int(unit[6])+3):
                    #    existingTE_inf[unit[5]]['start'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7],'1', '-')
                    #print unit[10], 'start', unit[6]
                    #for i in range(int(unit[7])-2, int(unit[7])+3):
                    #    existingTE_inf[unit[5]]['end'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '-')
                    #print unit[10], 'end', unit[7]
                    intact = 0
                    if int(unit[14]) == 1 and len(unit[12]) == 3:
                        unit[12] =re.sub(r'\(|\)', '', unit[12])
                        if int(unit[12]) == 0:
                            intact = 1
                    print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])), str(int(unit[7])), unit[11],unit[6],unit[7], intact, '-')
    ofile_RM.close()

#Chr3	not.give	RelocaTE_i	283493	283504	.	-	.	ID=repeat_Chr3_283493_283504;TSD=ATGCCATCAAGG;Note=Non-reference,
#not found in reference;Right_junction_reads:4;Left_junction_reads:1;Right_support_reads:4;Left_support_reads:5;
#Chr3	281479	284272	TE110112	+
def Overlap_TE_boundary(prefix, refte):
    data = defaultdict(str)
    final_gff   = '%s.gff' %(prefix)
    raw_gff     = '%s.raw.gff' %(prefix)
    clean_gff   = '%s.clean.gff' %(prefix)
    infile = '%s.overlap' %(prefix)
    outfile= '%s.remove.gff' %(prefix)
    os.system('bedtools window -w 10 -a %s -b %s > %s' %(final_gff, refte, infile))
    if not os.path.isfile(infile) or not os.path.getsize(infile) > 0:
        return 1 
    ofile  = open(outfile, 'w') 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                temp = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    if not attr == '':
                        attr = re.sub(r':', '=', attr)
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                if int(temp['Right_junction_reads']) == 0 or int(temp['Left_junction_reads']) == 0:
                    #support by one junction
                    #within 10 bp interval of intact TE boundary
                    #print >> ofile, '\t'.join(unit[:9])
                    if int(unit[3]) >= int(unit[10]) - 10 and int(unit[3]) <= int(unit[10]) + 10:
                        print >> ofile, '\t'.join(unit[:9])
                    elif int(unit[3]) >= int(unit[11]) - 10 and int(unit[3]) <= int(unit[11]) + 10:
                        print >> ofile, '\t'.join(unit[:9])
                    elif int(unit[4]) >= int(unit[10]) - 10 and int(unit[4]) <= int(unit[10]) + 10:
                        print >> ofile, '\t'.join(unit[:9])
                    elif int(unit[4]) >= int(unit[11]) - 10 and int(unit[4]) <= int(unit[11]) + 10:
                        print >> ofile, '\t'.join(unit[:9])
    ofile.close()
    if not os.path.isfile(outfile) or not os.path.getsize(outfile) > 0:
        return 1
    os.system('bedtools intersect -v -a %s -b %s > %s' %(final_gff, outfile, clean_gff))
    os.system('mv %s %s' %(final_gff, raw_gff))
    os.system('grep -v \"singleton\|insufficient_data\" %s > %s' %(clean_gff, final_gff))
    os.system('rm %s.overlap %s.remove.gff %s.clean.gff' %(prefix, prefix, prefix))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam')
    parser.add_argument('-t', '--te_fasta')
    parser.add_argument('-d', '--fq_dir')
    parser.add_argument('-g', '--genome_fasta')
    parser.add_argument('-r', '--reference_ins')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-s', '--size', default='500', type=int)
    parser.add_argument('-c', '--cpu', default='1', type=int)
    parser.add_argument('--len_cut_match', default='10', type=int)
    parser.add_argument('--len_cut_trim', default='10', type=int)
    parser.add_argument('--mismatch', default='0', type=int)
    parser.add_argument('--step', default='1234567', type=str)
    parser.add_argument('--run', dest='run', action='store_true')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    try:
        os.path.isfile(args.bam) and os.path.isfile(args.te_fasta) and os.path.isfile(args.genome_fasta)
    except:
        try:
            os.path.exists(args.fq_dir) and os.path.isfile(args.te_fasta) and os.path.isfile(args.genome_fasta)
        except:
            usage()
            sys.exit(2)

    #absolute path for all files, directory and scripts
    RelocaTE_bin = os.path.split(os.path.abspath(__file__))[0]
    reference    = os.path.abspath(args.genome_fasta)
    te_fasta     = os.path.abspath(args.te_fasta)
    mode         = 'bam'
    bam          = ''
    fastq_dir    = ''
    
    try: 
        os.path.isfile(args.bam)
        bam  = os.path.abspath(args.bam)
    except:
        try:
            if os.path.abspath(args.fq_dir):
                fastq_dir = os.path.abspath(args.fq_dir)
                mode      = 'fastq'
        except:
            usage()
            exit(2)

    print fastq_dir
    print mode

    #Prepare directory and script
    if args.outdir is None:
        args.outdir = '%s/RelocaTE_output' %(os.getcwd())
        createdir(args.outdir)
    else:
        args.outdir = os.path.abspath(args.outdir)
        createdir(args.outdir)

    #Default value of parameters
    if args.size is None:
        args.size = 500
    
    if args.cpu is None:
        args.cpu  = 1
    
    if args.step is None:
        args.step = '1234567'

    samtools = ''
    bedtools = ''
    bwa      = ''
    
    try:
        subprocess.check_output('which samtools', shell=True)
        samtools = subprocess.check_output('which samtools', shell=True)
        samtools = re.sub(r'\n', '', samtools)
    except:
        samtools = '/opt/samtools-0.1.16/samtools'

    try:
        subprocess.check_output('which bedtools', shell=True)
        bedtools = subprocess.check_output('which bedtools', shell=True)
        bedtools = re.sub(r'\n', '', bedtools)
    except:
        bedtools = '/opt/bedtools/2.17.0-25-g7b42b3b/bin//bedtools'

    try:
        subprocess.check_output('which bwa', shell=True)
        bwa = subprocess.check_output('which bwa', shell=True)
        bwa = re.sub(r'\n', '', bwa)
    except:
        bwa = '/opt/tyler/bin/bwa'
 
    try:
        subprocess.check_output('which blat', shell=True)
        blat = subprocess.check_output('which blat', shell=True)
        blat = re.sub(r'\n', '', blat)
    except:
        blat = '/usr/local/bin/blat'   
    
    try:
        subprocess.check_output('which seqtk', shell=True)
        seqtk = subprocess.check_output('which seqtk', shell=True)
        seqtk = re.sub(r'\n', '', seqtk)
    except:
        seqtk = '/rhome/cjinfeng/software/tools/seqtk-master//seqtk'

    #MSU_r7.fa.bwt
    if not os.path.isfile('%s.bwt' %(reference)):
        print 'Reference need to be indexed by bwa: %s' %(reference)
        exit()
 
    run_std = '%s/run.std' %(args.outdir)

    writefile('%s/regex.txt' %(args.outdir), '_1\t_2\t.unPaired\tUNK')
    createdir('%s/shellscripts' %(args.outdir))
    createdir('%s/repeat' %(args.outdir))
    createdir('%s/repeat/blat_output' %(args.outdir))
    createdir('%s/repeat/flanking_seq' %(args.outdir))
    createdir('%s/repeat/te_containing_fq' %(args.outdir))
    createdir('%s/repeat/te_only_read_portions_fa' %(args.outdir))

    shells = []
    #step0 existing TE blat
    shells_step0 = []
    reference_ins_flag = 'NONE'
    if args.reference_ins is None or args.reference_ins == '0':
        step0_file = '%s/shellscripts/step_0_do_not_call_reference_insertions' %(args.outdir)
        if '1' in list(args.step):
            writefile(step0_file, '')
    elif os.path.isfile(args.reference_ins):
        step0_file = '%s/shellscripts/step_0_te_annotation_provided' %(args.outdir)
        reference_ins_flag = args.reference_ins
        if '1' in list(args.step):
            writefile(step0_file, '')
    elif args.reference_ins == '1':
        createdir('%s/shellscripts/step_0' %(args.outdir))
        step0_file = '%s/shellscripts/step_0/step_0.existingTE_blat.sh' %(args.outdir)
        if '1' in list(args.step):
            shells.append('sh %s' %(step0_file))
            shells_step0.append('sh %s' %(step0_file))
            existingTE_blat = '%s %s %s %s/existingTE.blatout 1> %s/existingTE.blat.stdout' %(blat, reference, te_fasta, args.outdir, args.outdir)
            reference_ins_flag = '%s/existingTE.blatout' %(args.outdir)
            writefile(step0_file, existingTE_blat)

    #run job in this script
    if args.run and len(shells_step0) > 0 and '1' in list(args.step):
        single_run(shells_step0)

    #step1 format reference genome

    #step2 fastq to fasta
    shells_step2 = []
    fastas = defaultdict(lambda : str)
    if mode == 'fastq':
        fastqs = glob.glob('%s/*.f*q*' %(fastq_dir))
        step2_flag = 0
        step2_count= 0
        for fq in fastqs:
            fa    = ''
            if os.path.splitext(fq)[-1] == '.gz':
                fa    = '%s.fa' %(os.path.splitext(os.path.splitext(fq)[0])[0])
            else:
                fa    = '%s.fa' %(os.path.splitext(fq)[0])
            fastas[fa] = fq
            if not os.path.isfile(fa):
                createdir('%s/shellscripts/step_2' %(args.outdir))
                fq2fa = '%s seq -A %s > %s' %(seqtk, fq, fa)
                step2_file = '%s/shellscripts/step_2/%s.fq2fa.sh' %(args.outdir, step2_count)
                if '2' in list(args.step):
                    shells.append('sh %s' %(step2_file))
                    shells_step2.append('sh %s' %(step2_file))
                    writefile(step2_file, fq2fa)
                step2_flag == 1
                step2_count += 1
        if step2_flag == 0:
            step2_file = '%s/shellscripts/step_2_not_needed_fq_already_converted_2_fa' %(args.outdir)
            if '2' in list(args.step):
                writefile(step2_file, '')
        
        #run job in this script
        if args.run and len(shells_step2) > 0 and '2' in list(args.step):
            if int(args.cpu) == 1:
                single_run(shells_step2)
            else:
                mp_pool(shells_step2, int(args.cpu))
 
    elif mode == 'bam':
        #print 'Add module of obtaining reads from bam then prepare as fa files'
        cmd_step2 = []
        fastq_dir = '%s/repeat/fastq' %(args.outdir)
        createdir(fastq_dir)
        subbam = '%s/%s.subset.bam' %(fastq_dir, os.path.splitext(os.path.split(bam)[1])[0])
        fq1 = '%s/%s_1.fq' %(fastq_dir, os.path.splitext(os.path.split(bam)[1])[0])
        fq2 = '%s/%s_2.fq' %(fastq_dir, os.path.splitext(os.path.split(bam)[1])[0])
        fa1 = '%s.fa' %(os.path.splitext(fq1)[0])
        fa2 = '%s.fa' %(os.path.splitext(fq2)[0])
        fastas[fa1] = fq1
        fastas[fa2] = fq2
        if not os.path.isfile(subbam):
            cmd_step2.append('%s view -h %s | awk \'$5<60\' | samtools view -Shb - | samtools sort -m 500000000 -n - %s 2> %s' %(samtools, bam, os.path.splitext(subbam)[0], run_std))
        if not os.path.isfile(fq1) and not os.path.isfile(fq2):
            cmd_step2.append('%s bamtofastq -i %s -fq %s -fq2 %s 2> %s' %(bedtools, subbam, fq1, fq2, run_std))
        cmd_step2.append('%s seq -A %s > %s' %(seqtk, fq1, fa1))
        cmd_step2.append('%s seq -A %s > %s' %(seqtk, fq2, fa2))
        step2_flag = 0
        if not os.path.isfile(fa1) and not os.path.isfile(fa2):
            createdir('%s/shellscripts/step_2' %(args.outdir))
            step2_file = '%s/shellscripts/step_2/0.bam2fa.sh' %(args.outdir)
            if '2' in list(args.step):
                shells.append('sh %s' %(step2_file))
                shells_step2.append('sh %s' %(step2_file))
                writefile(step2_file, '\n'.join(cmd_step2))
        else:
            step2_file = '%s/shellscripts/step_2_not_needed_fq_already_converted_2_fa' %(args.outdir)
            if '2' in list(args.step):
                writefile(step2_file, '')

        #run job in this script
        if args.run and len(shells_step2) > 0 and '2' in list(args.step):
            single_run(shells_step2)   

    #step3 blat fasta to repeat
    shells_step3 = []
    step3_count = 0
    for fa in sorted(fastas.keys()):
        createdir('%s/shellscripts/step_3' %(args.outdir))
        fq      = fastas[fa]
        #fq      = '%s.fq' %(os.path.splitext(fa)[0]) if os.path.isfile('%s.fq' %(os.path.splitext(fa)[0])) else '%s.fastq' %(os.path.splitext(fa)[0])
        fa_prefix = os.path.split(os.path.splitext(fa)[0])[1]
        blatout = '%s/repeat/blat_output/%s.te_repeat.blatout' %(args.outdir, fa_prefix)
        blatstd = '%s/repeat/blat_output/blat.out' %(args.outdir)
        blatcmd = '%s -minScore=10 -tileSize=7 %s %s %s 1>>%s 2>>%s' %(blat ,te_fasta, fa, blatout, blatstd, blatstd)
        flank= '%s/repeat/flanking_seq/%s.te_repeat.flankingReads.fq' %(args.outdir, fa_prefix)
        trim = 'python %s/relocaTE_trim.py %s %s %s %s %s > %s' %(RelocaTE_bin, blatout, fq, args.len_cut_match, args.len_cut_trim, args.mismatch, flank)
        step3_file = '%s/shellscripts/step_3/%s.te_repeat.blat.sh' %(args.outdir, step3_count)
        if not os.path.isfile(blatout) or os.path.getsize(blatout) == 0:
            if '3' in list(args.step):
                shells.append('sh %s' %(step3_file))
                shells_step3.append('sh %s' %(step3_file))
                step3_cmds = '%s\n%s' %(blatcmd, trim)
                writefile(step3_file, step3_cmds)
                step3_count += 1
        elif not os.path.isfile(flank) or os.path.getsize(flank) == 0:
            if '3' in list(args.step):
                shells.append('sh %s' %(step3_file))
                shells_step3.append('sh %s' %(step3_file))
                step3_cmds = '%s' %(trim)
                writefile(step3_file, step3_cmds)
                step3_count += 1
       
    #run job in this script
    if args.run and len(shells_step3) > 0 and '3' in list(args.step):
        if int(args.cpu) == 1:
            single_run(shells_step3)
        else:
            mp_pool(shells_step3, int(args.cpu)) 


    #step4 align TE trimed reads to genome
    shells_step4 = []
    ref = os.path.split(os.path.splitext(reference)[0])[1]
    createdir('%s/shellscripts/step_4' %(args.outdir))
    step4_file= '%s/shellscripts/step_4/step_4.%s.repeat.align.sh' %(args.outdir, ref)
    if '4' in list(args.step):
        shells.append('sh %s' %(step4_file))
        shells_step4.append('sh %s' %(step4_file))
        step4_cmd = 'python %s/relocaTE_align.py %s %s/repeat %s %s %s/regex.txt repeat not.given 0' %(RelocaTE_bin, RelocaTE_bin, args.outdir, reference, fastq_dir, args.outdir)
        writefile(step4_file, step4_cmd)
    
    #run job in this script
    if args.run and len(shells_step4) > 0 and '4' in list(args.step):
        single_run(shells_step4)

   
    #existingTE bed file
    top_dir = '%s/repeat' %(os.path.abspath(args.outdir))
    #read existing TE from file
    r_te = re.compile(r'repeatmasker|rm|\.out', re.IGNORECASE)
    if os.path.isfile(reference_ins_flag) and os.path.getsize(reference_ins_flag) > 0:
        if r_te.search(reference_ins_flag):
            existingTE_RM_ALL(top_dir, reference_ins_flag)
    else:
        print 'Existing TE file does not exists or zero size'
 
    #step5 find insertions
    shells_step5 = []
    ids = fasta_id(reference)
    createdir('%s/shellscripts/step_5' %(args.outdir))
    step5_count = 0
    for chrs in ids:
        step5_cmd = 'python %s/relocaTE_insertionFinder.py %s/repeat/bwa_aln/%s.repeat.bwa.sorted.bam %s %s repeat %s/regex.txt not.give 100 %s 0 0 %s' %(RelocaTE_bin, args.outdir, ref, chrs, reference, args.outdir, reference_ins_flag, args.size)
        step5_file= '%s/shellscripts/step_5/%s.repeat.findSites.sh' %(args.outdir, step5_count)
        if '5' in list(args.step):
            shells.append('sh %s' %(step5_file))
            shells_step5.append('sh %s' %(step5_file))
            writefile(step5_file, step5_cmd)
            step5_count +=1
    
    #run job in this script
    if args.run and len(shells_step5) > 0 and '5' in list(args.step):
        if int(args.cpu) == 1:
            single_run(shells_step5)
        else:
            mp_pool(shells_step5, int(args.cpu))   


    #step6 find transposons on reference: reference only or shared
    shells_step6 = []
    createdir('%s/shellscripts/step_6' %(args.outdir))
    step6_count = 0
    if mode == 'fastq':
        for chrs in ids:
            step6_cmd = 'python %s/relocaTE_absenceFinder.py %s/repeat/bwa_aln/%s.repeat.bwa.sorted.bam %s %s repeat %s/regex.txt not.give 100 %s 0 0 %s' %(RelocaTE_bin, args.outdir, ref, chrs, reference, args.outdir, reference_ins_flag, args.size)
            step6_file= '%s/shellscripts/step_6/%s.repeat.absence.sh' %(args.outdir, step6_count)
            if '6' in list(args.step):
                shells.append('sh %s' %(step6_file))
                shells_step6.append('sh %s' %(step6_file))
                writefile(step6_file, step6_cmd)
                step6_count +=1
        
        #run job in this script
        if args.run and len(shells_step6) > 0 and '6' in list(args.step):
            if int(args.cpu) == 1:
                single_run(shells_step6)
            else:
                mp_pool(shells_step6, int(args.cpu))   
        
    elif mode == 'bam':
        pass
 
    #step7 characterize homozygous, heterozygous and somatic insertion
    shells_step7 = []
    createdir('%s/shellscripts/step_7' %(args.outdir))
    step7_cmd = []
    step7_cmd.append('cat %s/repeat/results/*.all_nonref_insert.gff > %s/repeat/results/ALL.all_nonref_insert.gff' %(args.outdir, args.outdir))
    step7_cmd.append('python %s/clean_false_positive.py --input %s/repeat/results/ALL.all_nonref_insert.gff --refte %s/existingTE.bed' %(RelocaTE_bin, args.outdir, top_dir))
    step7_cmd.append('cat %s/repeat/results/*.all_nonref_insert.txt | grep "^TE" -v > %s/repeat/results/ALL.all_nonref_insert.txt' %(args.outdir, args.outdir))
    step7_cmd.append('cat %s/repeat/results/*.all_ref_insert.txt > %s/repeat/results/ALL.all_ref_insert.txt' %(args.outdir, args.outdir))
    step7_cmd.append('cat %s/repeat/results/*.all_ref_insert.gff > %s/repeat/results/ALL.all_ref_insert.gff' %(args.outdir, args.outdir))
    step7_cmd.append('perl %s/characterizer.pl -s %s/repeat/results/ALL.all_nonref_insert.txt -b %s -g %s --samtools %s' %(RelocaTE_bin, args.outdir, bam, reference, samtools))
    step7_file= '%s/shellscripts/step_7/0.repeat.characterize.sh' %(args.outdir)
    if '7' in list(args.step):
        shells.append('sh %s' %(step7_file))
        shells_step7.append('sh %s' %(step7_file)) 
        writefile(step7_file, '\n'.join(step7_cmd))
   
    #run job in this script
    if args.run and len(shells_step7) > 0 and '7' in list(args.step):
        single_run(shells_step7)
    

    #write script, always write cmd to file
    writefile('%s/run_these_jobs.sh' %(args.outdir), '\n'.join(shells))

if __name__ == '__main__':
    main()

