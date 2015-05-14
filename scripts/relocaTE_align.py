#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import subprocess

def usage():
    test="name"
    message='''
python relocaTE_align.py scripts path genome regex_file TE exper bowtie2

Align flanking_fq to genome in unpaired and paired manner
    '''
    print message

def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def fastq2id(fastq, id_file):
    awk_cmd0 = 'cat %s | ' %(fastq)
    awk_cmd1 = 'awk \'{if(NR%4 == 1){print $1}}\' | '
    awk_cmd2 = 'sed \'s/^@//\' > %s' %(id_file)
    awk_cmd  = awk_cmd0 + awk_cmd1 + awk_cmd2
    #print awk_cmd
    os.system(awk_cmd)   

def readid(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data.append(unit[0])
    return data

#read_500_1005237/2:end:3
#read_500_1005237/2
def paired_id(fullread1_id, fullread2_id, fullreadu_id):
    #os.system('cp %s %s.temp' %(fullread1_id, fullread1_id))
    #os.system('cp %s %s.temp' %(fullread2_id, fullread2_id))
    #id1 = readid('%s.temp' %(fullread1_id))
    #id2 = readid('%s.temp' %(fullread2_id))
    id1 = readid(fullread1_id)
    id2 = readid(fullread2_id)
    idu = readid(fullreadu_id)
    ofile_id1 = open(fullread1_id, 'w')
    ofile_id2 = open(fullread2_id, 'w')
    r_id= re.compile(r'(.*)\:(start|end):\d+')
    for i in range(len(id1)):
        flag1 = 0
        flag2 = 0
        if r_id.search(id1[i]):
            id1[i] = r_id.search(id1[i]).groups(0)[0]
            flag1  = 1 
        if r_id.search(id2[i]):
            id2[i] = r_id.search(id2[i]).groups(0)[0]
            flag2  = 1
        if flag1 == 1 or flag2 == 1:
            print >> ofile_id1, id1[i]
            print >> ofile_id2, id2[i]
    for i in range(len(idu)):
        if r_id.search(idu[i]):
            idu[i] = r_id.search(idu[i]).groups(0)[0]
            if idu[i][-2:] == '/1':
                print >> ofile_id1, idu[i]
                print >> ofile_id2, '%s/2' %(idu[i][:-2])
            elif idu[i][-2:] == '/2':
                print >> ofile_id2, idu[i]
                print >> ofile_id1, '%s/1' %(idu[i][:-2])
            else:
                print >> ofile_id1, idu[i]
                print >> ofile_id2, idu[i]            

    ofile_id1.close()
    ofile_id2.close()

def paired_fq(fullread1_id, fullread2_id, fastq1, fastq2, fq1_0_temp, fq2_0_temp, fastq_dir, seqtk):
    fastqs = glob.glob('%s/*.f*q*' %(fastq_dir)) 
    suffix = ''
    s = re.compile(r'(\.f\w*?q.*?$)')
    if s.search(fastqs[0]):
        suffix = s.search(fastqs[0]).groups(0)[0]
    fq1_0  = '%s/%s' %(fastq_dir, re.sub(r'.te_repeat.flankingReads.fq', r'%s' %(suffix), os.path.split(fastq1)[1]))
    fq2_0  = '%s/%s' %(fastq_dir, re.sub(r'.te_repeat.flankingReads.fq', r'%s' %(suffix), os.path.split(fastq2)[1]))
    #fq1_0_temp = re.sub(r'.id$', r'.fq', fullread1_id)
    #fq2_0_temp = re.sub(r'.id$', r'.fq', fullread2_id)
    seqtk_cmd1 = '%s subseq %s %s > %s' %(seqtk, fq1_0, fullread1_id, fq1_0_temp)
    seqtk_cmd2 = '%s subseq %s %s > %s' %(seqtk, fq2_0, fullread2_id, fq2_0_temp)
    #print seqtk_cmd1
    #print seqtk_cmd2
    os.system(seqtk_cmd1)
    os.system(seqtk_cmd2)

def parse_regex(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\s+',line)
                unit[0] = re.sub(r'.(fq|fastq)','',unit[0])
                unit[1] = re.sub(r'.(fq|fastq)','',unit[1])
                unit[2] = re.sub(r'.(fq|fastq)','',unit[2])
                data = unit
    return data

def find_mate_pair_lib(path, mate_file):
    flanking_fq = defaultdict(lambda : defaultdict(lambda : str))
    mate_file_1          = mate_file[0]
    mate_file_2          = mate_file[1]
    mate_file_unpaired   = mate_file[2]
    files_1 = []
    files_2 = []
    files_unpaired = []
    s_1 = re.compile(r'%s\.' %(mate_file_1))
    s_2 = re.compile(r'%s\.' %(mate_file_2))
    s_u = re.compile(r'%s' %(mate_file_unpaired))
    #print '%s,%s,%s' %(mate_file_1, mate_file_2, mate_file_unpaired) 
    flanking_files = glob.glob('%s/flanking_seq/*flankingReads.fq' %(path))
    for file0 in flanking_files:
        if os.path.getsize(file0) == 0:
            continue
        print file0
        if s_u.search(file0):
            print 'unpaired'
            files_unpaired.append(file0)
        elif s_1.search(file0):
            print 'p1'
            files_1.append(file0)
        elif s_2.search(file0):
            print 'p2'
            files_2.append(file0)
    if len(files_1) >= 1 and len(files_2) >= 1:
        for i in range(len(files_1)):
            file1 = files_1[i]
            file1 = re.sub(r'%s\..*' %(mate_file_1), '', file1)
            for j in range(len(files_2)):
                file2 = files_2[j]
                file2 = re.sub(r'%s\..*' %(mate_file_2), '', file2)
                #print 'file1: %s; file2: %s' %(file1, file2)
                if file1 == file2:
                    flanking_fq[file1][1] = files_1[i]
                    flanking_fq[file1][2] = files_2[j]
                    print 'pair1: %s; pair2: %s' %(flanking_fq[file1][1], flanking_fq[file1][2])
                    if len(files_unpaired) >= 1:
                        for k in range(len(files_unpaired)):
                            fileu = files_unpaired[k]
                            fileu = re.sub(r'%s\..*' %(mate_file_unpaired), '', fileu)
                            flanking_fq[file1]['unpaired'] = files_unpaired[k]
    ## fastq files possibly not match to pattern _p1, _p2, .unPaired provided
    else:
        files_unpaired = glob.glob('%s/flanking_seq/*flankingReads.fq' %(path))
        for i in range(len(files_unpaired)):
            #print 'all unpaired: %s' %(files_unpaired[i])
            flanking_fq[files_unpaired[i]]['unpaired'] = files_unpaired[i]
    return flanking_fq


def bwa_run(path, genome_file, fastq, fq_name, target, readclass, bwa, samtools):
    cmd = []
    cmd.append('%s aln %s %s > %s/bwa_aln/%s.%s.bwa.%s.sai' %(bwa, genome_file, fastq, path, target, fq_name, readclass))
    cmd.append('%s samse %s %s/bwa_aln/%s.%s.bwa.%s.sai %s | %s view -bhS - > %s/bwa_aln/%s.%s.bwa.%s.bam' %(bwa, genome_file, path, target, fq_name, readclass, fastq, samtools, path, target, fq_name, readclass))
    cmd.append('rm %s/bwa_aln/*.bwa.*.sai' %(path))
    os.system('\n'.join(cmd))
    #print '\n'.join(cmd)

#view -Shb -F 4 - >
def bwa_run_paired(path, genome_file, fastq1, fastq2, fq_name, target, bwa, samtools):
    cmd = []
    cmd.append('%s aln %s %s > %s/bwa_aln/%s.%s.bwa.mates.sai' %(bwa, genome_file, fastq1, path, target, os.path.split(os.path.splitext(fastq1)[0])[1]))
    cmd.append('%s aln %s %s > %s/bwa_aln/%s.%s.bwa.mates.sai' %(bwa, genome_file, fastq2, path, target, os.path.split(os.path.splitext(fastq2)[0])[1]))
    cmd.append('%s sampe %s %s/bwa_aln/%s.%s.bwa.mates.sai %s/bwa_aln/%s.%s.bwa.mates.sai %s %s | %s view -bhS - > %s/bwa_aln/%s.%s.bwa.mates.bam' %(bwa, genome_file, path, target, os.path.splitext(os.path.split(fastq1)[1])[0], path, target, os.path.splitext(os.path.split(fastq2)[1])[0], fastq1, fastq2, samtools, path, target, fq_name))
    cmd.append('rm %s/bwa_aln/*.bwa.mates.sai' %(path))
    os.system('\n'.join(cmd))
    #print '\n'.join(cmd)

def map_reads_bwa(scripts, flanking_fq, path, genome_file, fastq_dir, target, bwa, samtools, seqtk):
    bwa_out_files = []
    bwa_out_files_f = []
    ##map reads with bowtie
    for file_pre in sorted(flanking_fq.keys()):
        #map reads as unpaired, treat all reads as single
        #for file_type in sorted(flanking_fq[file_pre].keys()):
        #    fastq   = flanking_fq[file_pre][file_type]
        #    fq_name = os.path.splitext(os.path.split(fastq)[1])[0]
        #    bwa_run(path, genome_file, fastq, fq_name, target, 'single')
        #    bwa_out_files.append('%s/bwa_aln/%s.%s.bwa.single.sam' %(path, target, fq_name))
        #map reads as paired, find paired and unpaired and map seperately
        print 'pre: %s' %(file_pre)
        if flanking_fq[file_pre].has_key(1) and flanking_fq[file_pre].has_key(2):
            fastq1  = flanking_fq[file_pre][1]
            fastq2  = flanking_fq[file_pre][2]
            print '%s\n%s' %(fastq1, fastq2)
            fq_name = os.path.splitext(os.path.split(fastq1)[1])[0]
            if int(os.path.getsize(fastq1)) > 0 and int(os.path.getsize(fastq2)) > 0:
                #prepare paired file fastq1.matched, fastq2.matched and *.unPaired.fq
                #need to rewrite this part, get all pairs that not matched to TE, which will be used as suporting reads
                #cmd = '%s/clean_pairs_memory.pl -1 %s -2 %s 1> %s/flanking_seq/%s.unPaired.fq 2>> %s/%s.stderr' 
                #%(scripts, fastq1, fastq2, path, fq_name, path, target)
                cmd = '%s/clean_pairs_memory.py --fq1 %s --fq2 %s --repeat %s/te_containing_fq --fq_dir %s --seqtk %s' %(scripts, fastq1, fastq2, path, fastq_dir, seqtk)
                print cmd
                os.system(cmd)
            match1 = '%s.matched' %(fastq1)
            match2 = '%s.matched' %(fastq2)
            unpaired = '%s/flanking_seq/%s.unPaired.fq' %(path, fq_name)
            if int(os.path.getsize(match1)) > 0 and int(os.path.getsize(match2)) > 0:
                #map paired-reads
                bwa_run_paired(path, genome_file, match1, match2, fq_name, target, bwa, samtools)
                bwa_out_files.append('%s/bwa_aln/%s.%s.bwa.mates.bam' %(path, target, fq_name))
            if int(os.path.getsize(unpaired) > 0):
                #map unpaired-reads
                bwa_run(path, genome_file, unpaired, fq_name, target, 'unPaired', bwa, samtools)
                bwa_out_files.append('%s/bwa_aln/%s.%s.bwa.unPaired.bam' %(path, target, fq_name))
            #get full reads of junction reads and their pairs
            #map these reads to genome and use perfect mapped reads as control for false junctions
            fullread1_id = '%s.fullreads.id' %(match1)
            fullread2_id = '%s.fullreads.id' %(match2)
            fullread1_fq = '%s.fullreads.fq' %(match1)
            fullread2_fq = '%s.fullreads.fq' %(match2)
            fullreadu_id = '%s.fullreads.id' %(unpaired)
            fastq2id(match1, fullread1_id)
            fastq2id(match2, fullread2_id)
            fastq2id(unpaired, fullreadu_id)
            paired_id(fullread1_id, fullread2_id, fullreadu_id)
            paired_fq(fullread1_id, fullread2_id, fastq1, fastq2, fullread1_fq, fullread2_fq, fastq_dir, seqtk)
            bwa_run_paired(path, genome_file, fullread1_fq, fullread2_fq, '%s.matched.fullreads' %(fq_name), target, bwa, samtools)
            bwa_out_files_f.append('%s/bwa_aln/%s.%s.bwa.mates.bam' %(path, target, '%s.matched.fullreads' %(fq_name))) 
        else:
        #paired not provided or not found, map reads as unpaired
            fastq   = flanking_fq[file_pre][file_type]
            fq_name = os.path.splitext(os.path.split(fastq)[1])[0]
            bwa_run(path, genome_file, fastq, fq_name, target, 'single', bwa, samtools)
            bwa_out_files.append('%s/bwa_aln/%s.%s.bwa.single.bam' %(path, target, fq_name))
   
    ##prepare merged files of bwa output
    #bam2merge = []
    #for i in range(len(bwa_out_files)):
    #    sam = bwa_out_files[i]
    #    bam = re.sub(r'sam', 'bam', sam)
    #    cmd = '%s view -h -bS %s > %s' %(samtools, sam, bam)
    #    os.system(cmd)
    #    bam2merge.append(bam)

    ##merge all bwa results into one file
    bam2merge  = bwa_out_files
    merged_bwa = '%s/bwa_aln/%s.repeat.bwa.bam' %(path, target)
    if len(bam2merge) > 1:
        cmd1  = '%s merge -f %s %s' %(samtools, merged_bwa, ' '.join(bam2merge))
        cmd2  = '%s sort %s %s.sorted' %(samtools, merged_bwa, os.path.splitext(merged_bwa)[0])
        cmd3  = '%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa)[0])
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
    elif len(bam2merge) == 1:
        os.system('cp %s %s' %(bam2merge[0], merged_bwa))
        os.system('%s sort %s %s.sorted' %(samtools, merged_bwa, os.path.splitext(merged_bwa)[0]))
        os.system('%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa)[0]))

    ##merge all bwa results of fullreads into one file
    bam2merge_f  = bwa_out_files_f
    merged_bwa_f = '%s/bwa_aln/%s.repeat.fullreads.bwa.bam' %(path, target)
    if len(bam2merge_f) > 1:
        cmd4  = '%s merge -f %s %s' %(samtools, merged_bwa_f, ' '.join(bam2merge_f))
        cmd5  = '%s sort %s %s.sorted' %(samtools, merged_bwa_f, os.path.splitext(merged_bwa_f)[0])
        cmd6  = '%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa_f)[0])
        os.system(cmd4)
        os.system(cmd5)
        os.system(cmd6)
    elif len(bam2merge_f) == 1:
        os.system('cp %s %s' %(bam2merge_f[0], merged_bwa_f))
        os.system('%s sort %s %s.sorted' %(samtools, merged_bwa_f, os.path.splitext(merged_bwa_f)[0]))
        os.system('%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa_f)[0]))


def bowtie_run(path, genome_file, fastq, fq_name, target, bowtie2, relax_align, bowtie_sam, readclass):
    if bowtie2 != 1 and bowtie_sam and relax_align != 1:
    #bowtie1 with sam output
        cmd1 = 'bowtie --sam --sam-nohead --sam-nosq -a -m 1 -v 3 -q %s.bowtie_build_index %s 1> %s/bowtie_aln/%s.%s.bowtie.%s.out 2>> %s/%s.stderr' %(genome_file, fastq, path, target, fq_name, readclass, path, target)
        os.system(cmd1)
    elif bowtie2 != 1 and bowtie_sam and relax_align == 1:
    ##bowtie1 with sam output, relax mapping
        cmd2 = 'bowtie --sam --sam-nohead --best -a -v 1 -q %s.bowtie_build_index %s 1> %s/bowtie_aln/%s.%s.bowtie.%s.out 2>> %s/%s.stderr' %(genome_file, fastq, path, target, fq_name, readclass, path, target)
        os.system(cmd2)
    elif bowtie2 == 1:
    ##bowtie2 -- need to get comparable -a -m1 -v3 arguments
        cmd3 = 'bowtie2 --sam-nohead --sam-nosq -x %s.bowtie_build_index -U %s 1> %s/bowtie_aln/%s.%s.bowtie.%s.out 2>> %s/%s.stderr' %(genome_file, fastq, path, target, fq_name, readclass, path, target)
        os.system(cmd3)
    else:
    #bowtie1 with bowtie output
        cmd4 = 'bowtie --best -q %s.bowtie_build_index %s 1> %s/bowtie_aln/%s.%s.bowtie.%s.out 2>> %s/%s.stderr' %(genome_file, fastq, path, target, fq_name, readclass, path, target)
        os.system(cmd4)
   
def bowtie_run_paired(path, genome_file, fastq1, fastq2, fq_name, target, bowtie2, relax_align, bowtie_sam):
    if bowtie2 != 1 and bowtie_sam and relax_align != 1:
    #bowtie1 with sam output
        cmd1 = 'bowtie --sam --sam-nohead --sam-nosq -a -m 1 -v 3 -q %s.bowtie_build_index -1 %s -2 %s 1> %s/bowtie_aln/%s.%s.bowtie.mates.out 2>> %s/%s.stderr' %(genome_file, fastq1, fastq2, path, target, fq_name, path, target)
        os.system(cmd1)
    elif bowtie2 != 1 and bowtie_sam and relax_align == 1:
    ##bowtie1 with sam output, relax mapping
        cmd2 = 'bowtie --sam --sam-nohead --best -a -v 1 -q %s.bowtie_build_index -1 %s -2 %s 1> %s/bowtie_aln/%s.%s.bowtie.mates.out 2>> %s/%s.stderr' %(genome_file, fastq1, fastq2, path, target, fq_name, path, target)
        os.system(cmd2)
    elif bowtie2 == 1:
    ##bowtie2 -- need to get comparable -a -m1 -v3 arguments
        cmd3 = 'bowtie2 --sam-nohead --sam-nosq -x %s.bowtie_build_index -1 %s -2 %s 1> %s/bowtie_aln/%s.%s.bowtie.mates.out 2>> %s/%s.stderr' %(genome_file, fastq1, fastq2, path, target, fq_name, path, target)
        os.system(cmd3)
    else:
    #bowtie1 with bowtie output
        cmd4 = 'bowtie --best -q %s.bowtie_build_index -1 %s -2 %s 1> %s/bowtie_aln/%s.%s.bowtie.mates.out 2>> %s/%s.stderr' %(genome_file, fastq, path, target, fq_name, path, target)
        os.system(cmd4)
 

def map_reads_bowtie(scripts, flanking_fq, path, genome_file, fastq_dir, target, bowtie2, relax_align, bowtie_sam):
    bowtie_out_files = []
    ##map reads with bowtie
    for file_pre in sorted(flanking_fq.keys()):
        #map reads as unpaired, treat all reads as single
        for file_type in sorted(flanking_fq[file_pre].keys()):
            fastq   = flanking_fq[file_pre][file_type]
            fq_name = os.path.splitext(os.path.split(fastq)[1])[0]
            bowtie_run(path, genome_file, fastq, fq_name, target, bowtie2, relax_align, bowtie_sam, 'single')
            bowtie_out_files.append('%s/bowtie_aln/%s.%s.bowtie.single.out' %(path, target, fq_name))
        #map reads as paired, find paired and unpaired and map seperately
        if flanking_fq[file_pre].has_key(1) and flanking_fq[file_pre].has_key(2):
            fastq1  = flanking_fq[file_pre][1]
            fastq2  = flanking_fq[file_pre][2]
            fq_name = os.path.splitext(os.path.split(fastq1)[1])[0]
            if int(os.path.getsize(fastq1)) > 0 and int(os.path.getsize(fastq2)) > 0:
                #prepare paired file fastq1.matched, fastq2.matched and *.unPaired.fq
                #need to rewrite this part, get all pairs that not matched to TE, which will be used as suporting reads
                #cmd = '%s/clean_pairs_memory.pl -1 %s -2 %s 1> %s/flanking_seq/%s.unPaired.fq 2>> %s/%s.stderr' %(scripts, fastq1, fastq2, path, fq_name, path, target)
                cmd = '%s/clean_pairs_memory.py --fq1 %s --fq2 %s --repeat %s/te_containing_fq --fq_dir %s' %(scripts, fastq1, fastq2, path, fastq_dir)
                os.system(cmd)
            match1 = '%s.matched' %(fastq1)
            match2 = '%s.matched' %(fastq2)
            unpaired = '%s/flanking_seq/%s.unPaired.fq' %(path, fq_name)
            if int(os.path.getsize(match1)) > 0 and int(os.path.getsize(match2)) > 0:
                #map paired-reads
                bowtie_run_paired(path, genome_file, match1, match2, fq_name, target, bowtie2, relax_align, bowtie_sam)
                bowtie_out_files.append('%s/bowtie_aln/%s.%s.bowtie.mates.out' %(path, target, fq_name))
                #map unpaired-reads
                bowtie_run(path, genome_file, unpaired, fq_name, target, bowtie2, relax_align, bowtie_sam, 'unPaired')
                bowtie_out_files.append('%s/bowtie_aln/%s.%s.bowtie.unPaired.out' %(path, target, fq_name))
   
    ##prepare merged files of bowtie output
    files2merge = ''
    filecount   = len(bowtie_out_files)
    if (int(filecount) > 50):
        ##too many files to merge
        big_files_2_merge = []
        count = 0
        for i in range(filecount)[::50]:
            count += 1
            files2merge = ' '.join(bowtie_out_files[i:i+50])
            cmd = 'cat %s > %s/bowtie_aln/%s.repeat.merged.bowtie.%s.temp' %(files2merge, path, target, count)
            os.system(cmd)
            big_files_2_merge.append('%s/bowtie_aln/%s.repeat.merged.bowtie.%s.temp' %(path, target, count))
        files2merge = ' '.join(big_files_2_merge)
    else:
        ##fewer files to merge
        files2merge = ' '.join(bowtie_out_files)

    ##merge all bowtie results into one file
    merged_bowtie = '%s/bowtie_aln/%s.repeat.bowtie.out' %(path, target)
    if files2merge != '':
        cmd = 'cat %s > %s' %(files2merge, merged_bowtie)
        os.system(cmd)

def main():
    if not len(sys.argv) == 9:
        usage()
        sys.exit(2)

    scripts    = sys.argv[1] #full path to scripts directory
    path       = sys.argv[2] #current/top/TE
    genome_file= sys.argv[3]
    fastq_dir  = sys.argv[4] #fastq_dir of raw reads, need to get the mates for TE-related reads
    regex_file = sys.argv[5] 
    TE         = sys.argv[6] #TE name, we use repeat for combined analysis
    exper      = sys.argv[7] #prefix for output
    bowtie2    = sys.argv[8] #use bowtie2 or not
    relax_align= 0  
    bowtie_sam = 1
    
    samtools = ''
    bwa      = ''
    seqtk    = ''    

    try:
        subprocess.check_output('which samtools', shell=True)
        samtools = subprocess.check_output('which samtools', shell=True)
        samtools = re.sub(r'\n', '', samtools)
    except:
        samtools = '/opt/samtools-0.1.16/samtools'

    try:
        subprocess.check_output('which bwa', shell=True)
        bwa = subprocess.check_output('which bwa', shell=True)
        bwa = re.sub(r'\n', '', bwa)
    except:
        bwa = '/opt/tyler/bin/bwa'
    
    
    try:
        subprocess.check_output('which seqtk', shell=True)
        seqtk = subprocess.check_output('which seqtk', shell=True)
        seqtk = re.sub(r'\n', '', seqtk)
    except:
        seqtk = '/rhome/cjinfeng/software/tools/seqtk-master//seqtk'

 
   
    #get the regelar expression patterns for mates and for the TE
    #when passed on the command line as an argument, even in single
    #quotes I lose special regex characters
    mate_file = parse_regex(regex_file)
    TSD                  = mate_file[3]
    #print TSD

    #get pairs of fastq and unpaired fastq files
    #Reads mapped to middle of TE need to removed and keep only their pairs if the pairs is not TE related
    flanking_fq = find_mate_pair_lib(path, mate_file)
    #for file_pre in sorted(flanking_fq.keys()):
    #    print 'file_pre keys:\n %s' %(file_pre)
    #    for k in sorted(flanking_fq[file_pre].keys()):
    #        print 'file:\n %s' %(flanking_fq[file_pre][k])

    #get directory path and file name in separate variables
    #map unpaired fastq and mate paired fastq
    target = os.path.splitext(os.path.split(genome_file)[1])[0]
    if 1:
        if not os.path.exists('%s/bwa_aln' %(path)):
            createdir('%s/bwa_aln' %(path))
        map_reads_bwa(scripts, flanking_fq, path, genome_file, fastq_dir, target, bwa, samtools, seqtk)
    else:
        if not os.path.exists('%s/bowtie_aln' %(path)):
            createdir('%s/bowtie_aln' %(path))
        map_reads_bowtie(scripts, flanking_fq, path, genome_file, fastq_dir, target, bowtie2, relax_align, bowtie_sam)

if __name__ == '__main__':
    main()

