#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import subprocess
import multiprocessing as mp
import pysam

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

##write content to new file
def writefile(outfile, lines):
    ofile = open(outfile, 'w')
    print >> ofile, lines
    ofile.close()


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data

def write_repeat_name_chr(read_chr_info, repeat_name, repeat_name_chr, read_order):
    read_order_pair = 1 if int(read_order) == 2 else 2
    ofile = open(repeat_name_chr, 'w')
    with open (repeat_name, 'r') as filehd: 
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                
                read_name = unit[0]
                if read_name[-2:] == '/1' or read_name[-2:] == '/2': read_name = read_name[:-2]
                if read_chr_info.has_key(read_name):
                    if read_chr_info[read_name].has_key(int(read_order)):
                        #by default set repeat to itself
                        print >> ofile, '%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], read_chr_info[read_name][int(read_order)])
                    elif read_chr_info[read_name].has_key(int(read_order_pair)):
                        #if self is not aligned, set to pair
                        print >> ofile, '%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], read_chr_info[read_name][int(read_order_pair)])
                    elif read_chr_info[read_name].has_key(0):
                        #set to unpaired
                        print >> ofile, '%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], read_chr_info[read_name][0])
                    else:
                        #should not happend, or very rare.
                        print >> ofile, '%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], 'NA')
                else:
                    #not in alignment, not usefull at all
                    #print >> ofile, '%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], 'NA')
                    pass
    ofile.close()

##read->read_order->chr
def get_read_chr(bam_file):
    read_chr = defaultdict(lambda : defaultdict(lambda : str()))
    r = re.compile(r'(.*):(start|end):(5|3)')
    fsam   = pysam.AlignmentFile(bam_file, 'rb')
    rnames = fsam.references
    rlens  = fsam.lengths
    for record in fsam.fetch(reference=None, until_eof = True):
        read_order = 0
        if record.is_read1:
            read_order = 1
        elif record.is_read2:
            read_order = 2
        qName    = record.query_name
        if r.search(qName): qName = r.search(qName).groups(0)[0]
        if qName[-2:] == '/1' or qName[-2:] == '/2': qName = qName[:-2]
        if not record.is_unmapped:
            #query inf
            #qName    = record.query_name
            #target inf
            tName    = rnames[record.reference_id]
            read_chr[qName][read_order] = tName
        else:
            #qName    = record.query_name
            read_chr[qName][read_order] = 'NA'
    return read_chr

def get_unpaired_info_chr(unpaired_bam, unpaired_info, unpaired_info_chr):
    info_dict = readtable(unpaired_info)
    ofile  = open(unpaired_info_chr, 'w')
    fsam   = pysam.AlignmentFile(unpaired_bam, 'rb')
    rnames = fsam.references
    rlens  = fsam.lengths
    for record in fsam.fetch(reference=None, until_eof = True):
        if not record.is_unmapped:
            #query inf
            qName    = record.query_name
            qLen     = int(record.query_length)
            qStart   = int(record.query_alignment_start)
            qEnd     = int(record.query_alignment_end) - 1
            #target inf
            tName    = rnames[record.reference_id]
            tLen     = int(rlens[record.reference_id])
            tStart   = int(record.reference_start)
            tEnd     = int(record.reference_end) - 1 
            if info_dict.has_key(qName):
                print >> ofile, '%s\t%s\t%s' %(qName, info_dict[qName], tName)
        else:
            qName    = record.query_name
            if info_dict.has_key(qName):
                print >> ofile, '%s\t%s\t%s' %(qName, info_dict[qName], 'NA')
    ofile.close()

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
    ##for paired-end junction reads, if one read is junction we output paired id/sequence to fullreads1 and fullreads2
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
    #for unpaired junction reads, if read is junction read we output paired-end id/sequence to fullread1 and fullread2
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
        if verbose > 3: print file0
        if s_u.search(file0):
            if verbose > 3: print 'unpaired'
            files_unpaired.append(file0)
        elif s_1.search(file0):
            if verbose > 3: print 'p1'
            files_1.append(file0)
        elif s_2.search(file0):
            if verbose > 3: print 'p2'
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
                    if verbose > 3: print 'pair1: %s; pair2: %s' %(flanking_fq[file1][1], flanking_fq[file1][2])
                    if len(files_unpaired) >= 1:
                        for k in range(len(files_unpaired)):
                            fileu = files_unpaired[k]
                            fileu = re.sub(r'%s\..*' %(mate_file_unpaired), '', fileu)
                            flanking_fq[file1]['unpaired'] = files_unpaired[k]
    ## fastq files possibly not match to pattern _p1, _p2, .unPaired provided
    else:
        files_unpaired = glob.glob('%s/flanking_seq/*flankingReads.fq' %(path))
        for i in range(len(files_unpaired)):
            if verbose > 3: print 'all unpaired: %s' %(files_unpaired[i])
            flanking_fq[files_unpaired[i]]['unpaired'] = files_unpaired[i]
    return flanking_fq


def bwa_run(path, genome_file, fastq, fq_name, target, readclass, bwa, samtools):
    cmd = []
    cmd.append('%s aln %s %s > %s/bwa_aln/%s.%s.bwa.%s.sai 2>>%s/bwa_aln/bwa.stderr' %(bwa, genome_file, fastq, path, target, fq_name, readclass, path))
    cmd.append('%s samse %s %s/bwa_aln/%s.%s.bwa.%s.sai %s | %s view -bhS - > %s/bwa_aln/%s.%s.bwa.%s.bam 2>>%s/bwa_aln/bwa.stderr' %(bwa, genome_file, path, target, fq_name, readclass, fastq, samtools, path, target, fq_name, readclass, path))
    cmd.append('rm %s/bwa_aln/%s.%s.bwa.%s.sai' %(path, target, fq_name, readclass))
    os.system('\n'.join(cmd))
    #print '\n'.join(cmd)

#view -Shb -F 4 - >
def bwa_run_paired(path, genome_file, fastq1, fastq2, fq_name, target, bwa, samtools):
    cmd = []
    cmd.append('%s aln %s %s > %s/bwa_aln/%s.%s.bwa.mates.sai 2>>%s/bwa_aln/bwa.stderr' %(bwa, genome_file, fastq1, path, target, os.path.split(os.path.splitext(fastq1)[0])[1], path))
    cmd.append('%s aln %s %s > %s/bwa_aln/%s.%s.bwa.mates.sai 2>>%s/bwa_aln/bwa.stderr' %(bwa, genome_file, fastq2, path, target, os.path.split(os.path.splitext(fastq2)[0])[1], path))
    cmd.append('%s sampe %s %s/bwa_aln/%s.%s.bwa.mates.sai %s/bwa_aln/%s.%s.bwa.mates.sai %s %s | %s view -bhS - > %s/bwa_aln/%s.%s.bwa.mates.bam 2>>%s/bwa_aln/bwa.stderr' %(bwa, genome_file, path, target, os.path.splitext(os.path.split(fastq1)[1])[0], path, target, os.path.splitext(os.path.split(fastq2)[1])[0], fastq1, fastq2, samtools, path, target, fq_name, path))
    cmd.append('rm %s/bwa_aln/%s.%s.bwa.mates.sai' %(path, target, os.path.split(os.path.splitext(fastq1)[0])[1]))
    cmd.append('rm %s/bwa_aln/%s.%s.bwa.mates.sai' %(path, target, os.path.split(os.path.splitext(fastq2)[0])[1]))
    os.system('\n'.join(cmd))
    #print '\n'.join(cmd)

def map_reads_bwa_mp_helper(args):
    return map_reads_bwa_mp_runner(*args)

#def map_reads_bwa_mp_runner(flanking_fq_list, scripts, path, genome_file, fastq_dir, target, bwa, samtools, seqtk, file_pre):
#def map_reads_bwa_mp_runner(seqtk, file_pre):
#    print 'CK: %s\t%s' %(seqtk, file_pre)

def map_reads_bwa_mp_runner(flanking_fq_list, scripts, path, genome_file, fastq_dir, target, bwa, samtools, seqtk, file_pre, mate_file):
    out_files   = []
    out_files_f = []
    if len(flanking_fq_list) == 2:
        fastq1  = flanking_fq_list[0]
        fastq2  = flanking_fq_list[1]
        if verbose > 3: print '%s\n%s' %(fastq1, fastq2)
        fq_name = os.path.splitext(os.path.split(fastq1)[1])[0]
        if int(os.path.getsize(fastq1)) > 0 and int(os.path.getsize(fastq2)) > 0:
            #prepare paired file fastq1.matched, fastq2.matched and *.unPaired.fq
            #need to rewrite this part, get all pairs that not matched to TE, which will be used as suporting reads
            #cmd = '%s/clean_pairs_memory.pl -1 %s -2 %s 1> %s/flanking_seq/%s.unPaired.fq 2>> %s/%s.stderr' 
            #%(scripts, fastq1, fastq2, path, fq_name, path, target)
            cmd = '%s/clean_pairs_memory.py --fq1 %s --fq2 %s --repeat %s/te_containing_fq --fq_dir %s --seqtk %s' %(scripts, fastq1, fastq2, path, fastq_dir, seqtk)
            if verbose > 3: print cmd
            os.system(cmd)
        match1 = '%s.matched' %(fastq1)
        match2 = '%s.matched' %(fastq2)
        unpaired = '%s/flanking_seq/%s.unPaired.fq' %(path, fq_name)
        unpaired_bam      = '%s/bwa_aln/%s.%s.bwa.unPaired.bam' %(path, target, fq_name)
        unpaired_info     = '%s/flanking_seq/%s.unPaired.info' %(path, fq_name)
        unpaired_info_chr = '%s/flanking_seq/%s.unPaired.info.chr' %(path, fq_name)
        mate_bam          = '%s/bwa_aln/%s.%s.bwa.mates.bam' %(path, target, fq_name)
        repeat_name_1     = '%s/te_containing_fq/%s%s.te_repeat.read_repeat_name.txt' %(path, os.path.split(file_pre)[1], mate_file[0])
        repeat_name_chr_1 = '%s/te_containing_fq/%s%s.te_repeat.read_repeat_name.chr.txt' %(path, os.path.split(file_pre)[1], mate_file[0])
        repeat_name_2     = '%s/te_containing_fq/%s%s.te_repeat.read_repeat_name.txt' %(path, os.path.split(file_pre)[1], mate_file[1])
        repeat_name_chr_2 = '%s/te_containing_fq/%s%s.te_repeat.read_repeat_name.chr.txt' %(path, os.path.split(file_pre)[1], mate_file[1])
        read_chr_info     = defaultdict(lambda : defaultdict(lambda : str()))
        if int(os.path.getsize(match1)) > 0 and int(os.path.getsize(match2)) > 0:
            #map paired-reads
            bwa_run_paired(path, genome_file, match1, match2, fq_name, target, bwa, samtools)
            out_files.append('%s/bwa_aln/%s.%s.bwa.mates.bam' %(path, target, fq_name))
            read_chr_info.update(get_read_chr(mate_bam))
        if int(os.path.getsize(unpaired) > 0):
            #map unpaired-reads
            bwa_run(path, genome_file, unpaired, fq_name, target, 'unPaired', bwa, samtools)
            get_unpaired_info_chr(unpaired_bam, unpaired_info, unpaired_info_chr)
            out_files.append('%s/bwa_aln/%s.%s.bwa.unPaired.bam' %(path, target, fq_name))
            read_chr_info.update(get_read_chr(unpaired_bam))
        write_repeat_name_chr(read_chr_info, repeat_name_1, repeat_name_chr_1, 1)
        write_repeat_name_chr(read_chr_info, repeat_name_2, repeat_name_chr_2, 2)
        read_chr_info.clear()
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
        out_files_f.append('%s/bwa_aln/%s.%s.bwa.mates.bam' %(path, target, '%s.matched.fullreads' %(fq_name))) 
    else:
        #paired not provided or not found, map reads as unpaired
        fastq   = flanking_fq_list[0]
        fq_name = os.path.splitext(os.path.split(fastq)[1])[0]
        bwa_run(path, genome_file, fastq, fq_name, target, 'single', bwa, samtools)
        out_files.append('%s/bwa_aln/%s.%s.bwa.single.bam' %(path, target, fq_name))
    
    return [out_files, out_files_f]


##run function with parameters using multiprocess of #cpu
def multiprocess_pool(parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(map_reads_bwa_mp_helper, tuple(parameters))
    collect_list = []
    for x in imap_it:
        #print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list
    #return 1

def map_reads_bwa(scripts, flanking_fq, path, genome_file, fastq_dir, target, bwa, samtools, seqtk, cpu, mate_file):
    bwa_out_files = []
    bwa_out_files_f = []
    ##map reads with bwa
    #hg18.p00.chr1_22_reads_10X_100_500_1.te_repeat.flankingReads.bwa.mates.bam
    parameters = []
    test_bam = '%s/bwa_aln/%s.%s%s.te_repeat.flankingReads.bwa.mates.bam' %(path, target, os.path.split(flanking_fq.keys()[0])[1], mate_file[0])
    if verbose > 0: print 'testing if bam exists: %s' %(test_bam)
        
    for file_pre in sorted(flanking_fq.keys()):
        #map reads as unpaired, treat all reads as single
        #for file_type in sorted(flanking_fq[file_pre].keys()):
        #    fastq   = flanking_fq[file_pre][file_type]
        #    fq_name = os.path.splitext(os.path.split(fastq)[1])[0]
        #    bwa_run(path, genome_file, fastq, fq_name, target, 'single')
        #    bwa_out_files.append('%s/bwa_aln/%s.%s.bwa.single.sam' %(path, target, fq_name))
        #map reads as paired, find paired and unpaired and map seperately
        if verbose > 3: print 'pre: %s' %(file_pre)
        flanking_fq_list = []
        if flanking_fq[file_pre].has_key(1) and flanking_fq[file_pre].has_key(2):
            flanking_fq_list = [flanking_fq[file_pre][1], flanking_fq[file_pre][2]]
        else:
            flanking_fq_list = [flanking_fq[file_pre]['unpaired']]
        parameters.append([flanking_fq_list, scripts, path, genome_file, fastq_dir, target, bwa, samtools, seqtk, file_pre, mate_file])
        #parameters.append([seqtk, file_pre])

    #for pm in parameters:
    #    print 'flanking_fq_list:', pm[0]

    ##mp runner
    if not os.path.isfile(test_bam):
        if verbose > 0: print 'bam not exists, preceed with bwa to map the reads'
        collect_files = multiprocess_pool(parameters, cpu)
        for run_files_list in collect_files:
            bwa_out_files.extend(run_files_list[0]) 
            bwa_out_files_f.extend(run_files_list[1])
    else:
        if verbose > 0: print 'bam exists, merge and sort bam files'
        bwa_out_files   = glob.glob('%s/bwa_aln/*.te_repeat.flankingReads.bwa.*.bam' %(path))
        bwa_out_files_f = glob.glob('%s/bwa_aln/*.te_repeat.flankingReads.matched.fullreads.bwa.*.bam' %(path))

    ##merge all bwa results into one file
    bam2merge  = bwa_out_files
    merged_bwa = '%s/bwa_aln/%s.repeat.bwa.bam' %(path, target)
    print 'mergeing bam file: %s/%s files' %(len(bam2merge), len(bwa_out_files))
    if len(bam2merge) > 1:
        cmd1  = '%s merge -f %s %s' %(samtools, merged_bwa, ' '.join(bam2merge))
        cmd2  = '%s sort %s %s.sorted' %(samtools, merged_bwa, os.path.splitext(merged_bwa)[0])
        cmd3  = '%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa)[0])
        #print '%s\n%s\n%s' %(cmd1, cmd2, cmd3)
        cmd_sh = '%s.sh' %(merged_bwa)
        writefile(cmd_sh, '%s\n%s\n%s' %(cmd1, cmd2, cmd3))
        os.system('bash %s' %(cmd_sh))
        #os.system(cmd1)
        #os.system(cmd2)
        #os.system(cmd3)
    elif len(bam2merge) == 1:
        os.system('cp %s %s' %(bam2merge[0], merged_bwa))
        os.system('%s sort %s %s.sorted' %(samtools, merged_bwa, os.path.splitext(merged_bwa)[0]))
        os.system('%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa)[0]))

    ##merge all bwa results of fullreads into one file
    bam2merge_f  = bwa_out_files_f
    merged_bwa_f = '%s/bwa_aln/%s.repeat.fullreads.bwa.bam' %(path, target)
    print 'mergeing fullread bam file: %s/%s files' %(len(bam2merge_f), len(bwa_out_files_f))
    if len(bam2merge_f) > 1:
        cmd4  = '%s merge -f %s %s' %(samtools, merged_bwa_f, ' '.join(bam2merge_f))
        cmd5  = '%s sort %s %s.sorted' %(samtools, merged_bwa_f, os.path.splitext(merged_bwa_f)[0])
        cmd6  = '%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa_f)[0])
        #print '%s\n%s\n%s' %(cmd4, cmd5, cmd6)
        cmd_sh = '%s.sh' %(merged_bwa_f)
        writefile(cmd_sh, '%s\n%s\n%s' %(cmd4, cmd5, cmd6))
        os.system('bash %s' %(cmd_sh))
        #os.system(cmd4)
        #os.system(cmd5)
        #os.system(cmd6)
    elif len(bam2merge_f) == 1:
        os.system('cp %s %s' %(bam2merge_f[0], merged_bwa_f))
        os.system('%s sort %s %s.sorted' %(samtools, merged_bwa_f, os.path.splitext(merged_bwa_f)[0]))
        os.system('%s index %s.sorted.bam' %(samtools, os.path.splitext(merged_bwa_f)[0]))

    ##merge info_chr file and split into chr based files
    info_chr_fh    = {}
    info_chr_files = glob.glob('%s/flanking_seq/*.unPaired.info.chr' %(path))
    for info_chr_file in info_chr_files:
        with open (info_chr_file, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2: 
                    unit = re.split(r'\t',line)
                    if not info_chr_fh.has_key(unit[2]):
                        info_chr_fh[unit[2]] = open('%s/flanking_seq/%s.flankingReads.unPaired.info' %(path, unit[2]), 'w')
                        print >> info_chr_fh[unit[2]], line
                    else:
                        print >> info_chr_fh[unit[2]], line
    for temp_chr in info_chr_fh:
        info_chr_fh[temp_chr].close()

    ##merge repeat_name_chr file and split into chr base files
    repeat_chr_fh    = {}
    repeat_chr_files = glob.glob('%s/te_containing_fq/*.te_repeat.read_repeat_name.chr.txt' %(path))
    for repeat_chr_file in repeat_chr_files:
        with open (repeat_chr_file, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2: 
                    unit = re.split(r'\t',line)
                    if not repeat_chr_fh.has_key(unit[3]):
                        repeat_chr_fh[unit[3]] = open('%s/te_containing_fq/%s.read_repeat_name.split.txt' %(path, unit[3]), 'w')
                        print >> repeat_chr_fh[unit[3]], line
                    else:
                        print >> repeat_chr_fh[unit[3]], line
    for temp_chr in repeat_chr_fh:
        repeat_chr_fh[temp_chr].close() 

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
    if not len(sys.argv) == 14:
        usage()
        sys.exit(2)

    global verbose
    scripts    = sys.argv[1] #full path to scripts directory
    path       = sys.argv[2] #current/top/TE
    genome_file= sys.argv[3]
    fastq_dir  = sys.argv[4] #fastq_dir of raw reads, need to get the mates for TE-related reads
    regex_file = sys.argv[5] 
    TE         = sys.argv[6] #TE name, we use repeat for combined analysis
    exper      = sys.argv[7] #prefix for output
    bowtie2    = sys.argv[8] #use bowtie2 or not
    cpu        = sys.argv[9]
    verbose    = int(sys.argv[10]) #0, 1, 2, 3, 4 
    relax_align= 0  
    bowtie_sam = 1
    
    samtools = sys.argv[11]
    bwa      = sys.argv[12]
    seqtk    = sys.argv[13]
   
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
        map_reads_bwa(scripts, flanking_fq, path, genome_file, fastq_dir, target, bwa, samtools, seqtk, cpu, mate_file)
    else:
        if not os.path.exists('%s/bowtie_aln' %(path)):
            createdir('%s/bowtie_aln' %(path))
        map_reads_bowtie(scripts, flanking_fq, path, genome_file, fastq_dir, target, bowtie2, relax_align, bowtie_sam)

if __name__ == '__main__':
    main()

