#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import gzip
import pysam

def usage():
    test="name"
    message='''
python relocaTE_trim.py te_repeat.blatout 500_1.fq 10 10 0

sys.argv[1]: blatout or bam
sys.argv[2]: fastq
sys.argv[3]: length cutoff for matched region
sys.argv[4]: length cutoff for trimmed fastq
sys.argv[5]: mismatch allowance: 0,1,2,3 mismatch allowed on reads and repeat alignment

Parse blat or bam file. Write TE matched reads and their pairs into files.
*.te_repeat.flankingReads.fq: trimmed reads without TE and reads matched to the middle of TE
*.ContainingReads.fq: all reads have matches to TE
*.five_prime.fa/*.three_prime.fa: TE proportion of reads that matched to TE
*.potential_tandemInserts_containing_reads.list.txt: reads have two or more matches on TE and on same strand. Tandem insertions and Highly false positive.
    '''
    print message

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    for i in range(len(bases)):
        bases[i] = complement[bases[i]] if complement.has_key(bases[i]) else bases[i]
    return ''.join(bases)

def reverse_complement(seq):
    return complement(seq[::-1])


def parse_align_blat(infile, tandem, verbose):
    coord = defaultdict(lambda : defaultdict(lambda : str))
    ##align_file
    ofile = open(tandem, 'w')
    with open (infile, 'r') as filehd:
        for i in range(5):
            next(filehd)
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            match    = int(unit[0])
            mismatch = int(unit[1])
            strand   = unit[8]
            qName    = unit[9]
            qLen     = int(unit[10])
            qStart   = int(unit[11])
            qEnd     = int(unit[12]) - 1
            tName    = unit[13]
            tLen     = int(unit[14])
            tStart   = int(unit[15])
            #get all values into 1st base = 0 postion notation
            tEnd     = int(unit[16]) - 1 

            #filter bad alignment with two many or large gaps
            #int(unit[17]) >= 3, fewer than 3 matching blocks
            if int(unit[4]) > 1 or int(unit[5]) > 3 or int(unit[6]) > 1 or int(unit[7]) > 3 or int(unit[17]) >= 3:
                continue
                
            #boundary = 1 if int(qStart) == 0 or int(qEnd) + 1 == int(qLen) else 0
            addRecord = 0
            boundary_qry_left = 0
            boundary_tar_left = 0
            boundary_qry_right = 0
            boundary_tar_right = 0
            if int(qStart) == 0 or int(qStart) <= 2:
                boundary_qry_left  = 1
            if int(qEnd) + 1 == int(qLen) or int(qEnd) >= int(qLen) - 3:
                boundary_qry_right = 1 
            if int(tStart) == 0 or int(tStart) <= 2:
                boundary_tar_left  = 1
            if int(tEnd) + 1 == int(tLen) or int(tEnd) >= int(tLen) - 3:
                boundary_tar_right = 1
            #max boundary should be 2: 1. match one read end and one repeat end; 2. match two read end and internal of repeat
            #we expect more boundary and compare match length when having equal number of boundary
            boundary = boundary_qry_left + boundary_tar_left + boundary_qry_right + boundary_tar_right
            
            if verbose >= 3: print >> sys.stderr, qName, qLen, qStart, qEnd, tName, tLen, tStart, tEnd, match, mismatch, boundary
            if coord.has_key(qName):
                ##keep the best match to TE
                if int(boundary) > int(coord[qName]['boundary']):
                        addRecord = 1
                elif int(boundary) == int(coord[qName]['boundary']):
                    if int(match) > int(coord[qName]['match']):
                        addRecord = 1
                    else:
                        addRecord = 0
                else:
                    addRecord = 0
            else:
                addRecord = 1 
            #filter bad alignment with two many or large gaps
            #if int(unit[4]) > 1 or int(unit[5]) > 3 or int(unit[6]) > 1 or int(unit[7]) > 3:
            #    addRecord = 0
            #print >> sys.stderr, qName, qStart, qEnd, match, addRecord
            if addRecord == 1:
                coord[qName]['match']    = match
                coord[qName]['len']      = qLen
                coord[qName]['start']    = qStart
                coord[qName]['end']      = qEnd
                coord[qName]['tLen']     = tLen
                coord[qName]['mismatch'] = mismatch
                coord[qName]['strand']   = strand
                coord[qName]['tName']    = tName
                coord[qName]['tStart']   = tStart
                coord[qName]['tEnd']     = tEnd
                coord[qName]['boundary'] = boundary
                if verbose >= 3: print >> sys.stderr, qName, qLen, qStart, qEnd, tName, tLen, tStart, tEnd, match, mismatch, boundary
    ofile.close()
    return coord

def convert_tag(tag):
    tags = {}
    for t in tag:
        tags[t[0]] = t[1]
    return tags

def update_coord(header1, header, coord):
    coord[header]['start']    = int(coord[header1]['start'])
    coord[header]['len']      = int(coord[header1]['len'])
    coord[header]['end']      = int(coord[header1]['end'])
    coord[header]['tName']    = coord[header1]['tName']
    coord[header]['tStart']   = int(coord[header1]['tStart'])
    coord[header]['tEnd']     = int(coord[header1]['tEnd'])
    coord[header]['tLen']     = int(coord[header1]['tLen'])
    coord[header]['mismatch'] = int(coord[header1]['mismatch'])
    coord[header]['match']    = int(coord[header1]['match'])
    coord[header]['strand']   = coord[header1]['strand']
    coord[header]['boundary'] = int(coord[header1]['boundary'])
 
    del coord[header1]['start']
    del coord[header1]['len']
    del coord[header1]['end']
    del coord[header1]['tName']
    del coord[header1]['tStart'] 
    del coord[header1]['tEnd']
    del coord[header1]['tLen']
    del coord[header1]['mismatch']
    del coord[header1]['match']
    del coord[header1]['strand']
    del coord[header1]['boundary']

def parse_align_bwa(infile, tandem, verbose):
    coord = defaultdict(lambda : defaultdict(lambda : str))
    ##align_file
    ofile = open(tandem, 'w')
    fsam = pysam.AlignmentFile(infile, 'rb')
    rnames = fsam.references
    rlens  = fsam.lengths
    qlen_c = 0
    for record in fsam.fetch(reference=None, until_eof = True):
        if not record.is_unmapped:
            #print >> sys.stderr, record.query_length, record.query_alignment_start, record.query_alignment_end
            #print >> sys.stderr, record.reference_start, record.reference_end
            #query inf
            qName    = record.query_name
            qLen     = int(record.query_length)
            qStart   = int(record.query_alignment_start)
            qEnd     = 0
            try:
                qEnd     = int(record.query_alignment_end) - 1
            except:
                continue
            if qLen == 0:
                qLen = qlen_c
                qEnd = qLen + qEnd
            else:
                qlen_c = qLen
            #target inf
            tName    = rnames[record.reference_id]
            tLen     = int(rlens[record.reference_id])
            tStart   = int(record.reference_start)
            tEnd     = int(record.reference_end) - 1
            #match and mismatch
            tag      = record.tags if record.tags else []
            tags     = convert_tag(tag)
            #mismatch = int(tags['NM'])
            #match    = int(record.query_alignment_length) - mismatch
            match    = 0
            ins0     = 0
            del0     = 0
            for (key, length) in record.cigartuples:
                #print key, length
                if int(key) == 0:
                    match += length
                elif int(key) == 1:
                    ins0 += length
                elif int(key) == 2:
                    del0 += length
            mismatch= int(tags['NM']) - int(ins0) - int(del0)
            match   = match - mismatch
            #strand, flag is 0 is read if read is unpaired and mapped to plus strand
            strand   = ''
            flag     = record.flag
            if int(flag) == 0:
                strand = '+'
            else:
                strand = '-' if record.is_reverse else '+'
            #update data
            #boundary = 1 if int(qStart) == 0 or int(qEnd) + 1 == int(qLen) else 0 
            addRecord = 0
            boundary  = 0
            boundary_qry_left = 0
            boundary_tar_left = 0
            boundary_qry_right = 0
            boundary_tar_right = 0
            if int(qStart) == 0 or int(qStart) <= 2:
                boundary_qry_left  = 1
            if int(qEnd) + 1 == int(qLen) or int(qEnd) >= int(qLen) - 3:
                boundary_qry_right = 1 
            if int(tStart) == 0 or int(tStart) <= 2:
                boundary_tar_left  = 1
            if int(tEnd) + 1 == int(tLen) or int(tEnd) >= int(tLen) - 3:
                boundary_tar_right = 1
            #max boundary should be 2: 1. match one read end and one repeat end; 2. match two read end and internal of repeat
            #we expect more boundary and compare match length when having equal number of boundary
            boundary = boundary_qry_left + boundary_tar_left + boundary_qry_right + boundary_tar_right
            if verbose >= 3: print >> sys.stderr, qName, qLen, qStart, qEnd, tName, tLen, tStart, tEnd, match, mismatch, boundary
            if coord.has_key(qName):
                ##keep the best match to TE
                if int(boundary) > int(coord[qName]['boundary']):
                        addRecord = 1
                elif int(boundary) == int(coord[qName]['boundary']):
                    if int(match) > int(coord[qName]['match']):
                        addRecord = 1
                    else:
                        addRecord = 0
                else:
                    addRecord = 0
            else:
                addRecord = 1
             
            #final data
            #print >> sys.stderr, qName, qStart, qEnd, match, addRecord
            if addRecord == 1:
                coord[qName]['match']    = match
                coord[qName]['len']      = qLen
                coord[qName]['start']    = qStart
                coord[qName]['end']      = qEnd
                coord[qName]['tLen']     = tLen
                coord[qName]['mismatch'] = mismatch
                coord[qName]['strand']   = strand
                coord[qName]['tName']    = tName
                coord[qName]['tStart']   = tStart
                coord[qName]['tEnd']     = tEnd
                coord[qName]['boundary'] = boundary
                if verbose >= 3: print >> sys.stderr, qName, qLen, qStart, qEnd, tName, tLen, tStart, tEnd, match, mismatch, boundary
    ofile.close()
    return coord


def main():
    if not len(sys.argv) == 6:
        usage()
        sys.exit(2)

    align_file   = sys.argv[1]
    fq_file1     = sys.argv[2]
    len_cutoff_m = int(sys.argv[3])
    len_cutoff_l = int(sys.argv[4])
    mismatch_allowance = int(sys.argv[5])
    verbose  = 0
    align_type = 'bam' if os.path.splitext(align_file)[1].replace(r'.', '') == 'bam' else 'blat'
  
    #set output directories and files
    align_path = os.path.abspath(os.path.split(align_file)[0])
    file_name  = os.path.split(align_file)[1]
    te_path    = re.split(r'/', align_path)
    del te_path[-1]
    out_fq_path= '%s/te_containing_fq' %('/'.join(te_path))
    out_fa_path= '%s/te_only_read_portions_fa' %('/'.join(te_path))
    tandem_file= '%s/%s.potential_tandemInserts_containing_reads.list.txt' %(out_fq_path, os.path.splitext(file_name)[0])
    read_repeat_file = '%s/%s.read_repeat_name.txt' %(out_fq_path, os.path.splitext(file_name)[0]) 
    #parse align
    TE = 'unspecified'
    FA = 'unspecified'
    coord = defaultdict(lambda : defaultdict(lambda : str))
    if align_type == 'blat':
        coord = parse_align_blat(align_file, tandem_file, verbose)
        s = re.compile(r'(\S+)\.te_(\S+)\.blat*')
        m = s.search(file_name)
        if m:
            FA = m.groups(0)[0]
            TE = m.groups(0)[1]
    elif align_type == 'bam':
        coord = parse_align_bwa(align_file, tandem_file, verbose)   
        s = re.compile(r'(\S+)\.te_(\S+)\.bam')
        m = s.search(file_name)
        if m:
            FA = m.groups(0)[0]
            TE = m.groups(0)[1] 

    #outfiles
    #TE = 'unspecified'
    #FA = 'unspecified' 
    #s = re.compile(r'(\S+)\.te_(\S+)\.blatout')
    #m = s.search(file_name)
    #if m:
    #    FA = m.groups(0)[0]
    #    TE = m.groups(0)[1]
    outfq  = 0
    outte5 = 0
    outte3 = 0
    
    if os.path.isfile('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE)) and int(os.path.getsize('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE))) > 0:
        outfq = 1        
    if os.path.isfile('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE)) and int(os.path.getsize('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE))) > 0:
        outte5 = 1
    if os.path.isfile('%s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE)) and int(os.path.getsize('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE))) > 0:
        outte3 = 1
    #print 'check1: %s\t%s\t%s' %(str(outfq), str(outte5), str(outte3))
    if outfq + outte5 + outte3 == 3:
        exit(0)
    else:
        if os.path.isfile('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE)):
            cmd1 = 'rm %s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE)
            os.system(cmd1)
        if os.path.isfile('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE)):
            cmd2 = 'rm %s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE)
            os.system(cmd2)
        if os.path.isfile('%s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE)):
            cmd3 = 'rm %s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE)
            os.system(cmd3)
        ofile_fq = open('%s/%s.te_%s.ContainingReads.fq' %(out_fq_path, FA, TE), 'w')
        ofile_te5= open('%s/%s.te_%s.five_prime.fa' %(out_fa_path, FA, TE), 'w')
        ofile_te3= open('%s/%s.te_%s.three_prime.fa' %(out_fa_path, FA, TE), 'w')
        ofile_rr = open(read_repeat_file, 'w')
        ##go thru each fq record in the fq files. if the name of the seq is in the blat file
        ##trim the seq
        #f_type = 'rb' if os.path.splitext(fq_file1)[-1] == '.gz' else 'r'
        filehd = ''
        if os.path.splitext(fq_file1)[-1] == '.gz':
            filehd = gzip.open (fq_file1, 'rb')
        else:
            filehd = open (fq_file1, 'r')
        if 1:
            for line in filehd:
                line = line.rstrip()
                header = line[1:]
                header = re.split(r'\s+', header)[0]
                rl_name= header
                seq    = filehd.next().rstrip()
                qualh  = filehd.next().rstrip()
                qual   = filehd.next().rstrip()
                #bam only
                header1 = ''
                if align_type == 'bam' and (header[-2:] == '/1' or header[-2:] == '/2'):
                    header1  = header[:-2]
                    if coord.has_key(header1):
                        update_coord(header1, header, coord)
 
                if coord.has_key(header):
                    start    = int(coord[header]['start'])
                    length   = int(coord[header]['len'])
                    end      = int(coord[header]['end'])
                    tName    = coord[header]['tName']
                    tStart   = int(coord[header]['tStart'])
                    tEnd     = int(coord[header]['tEnd'])
                    tLen     = int(coord[header]['tLen'])
                    mismatch = int(coord[header]['mismatch'])
                    match    = int(coord[header]['match'])
                    strand   = coord[header]['strand']
                   
                    #want to cut and keep anything not matching to database TE
                    trimmed_seq  = ''
                    trimmed_qual = ''
                    if verbose >= 3: print >> sys.stderr, 'Check:', header, start, end, length, tName, tStart, tEnd, tLen, mismatch, match, strand
                    #print 'check2: %s\t%s\t%s'  %(str(tStart), str((length - (match + mismatch))), str((mismatch/(match + mismatch))))
                    ##query read overlaps 5' end of database TE & trimmed seq > cutoff
                    #int(start) <= 2 or int(end) >= int(length) - 3, we need the reads mapped boundary to align with te
                    #if tStart == 0 and (int(start) <= 2 or int(end) >= int(length) - 3) and (length - (match + mismatch)) > len_cutoff and (float(mismatch)/(float(match) + float(mismatch))) <= mismatch_allowance:
                    #if tStart <= 2 and (int(start) <= 2 or int(end) >= int(length) - 3) and (length - (match + mismatch)) >= len_cutoff_m and int(mismatch)  <= int(mismatch_allowance):
                    if tStart <= 2 and (int(start) <= 2 or int(end) >= int(length) - 3) and (match + mismatch) >= len_cutoff_m and (length - (match + mismatch)) >= len_cutoff_l and int(mismatch)  <= int(mismatch_allowance):
                        tS = int(tStart) + 1
                        tE = int(tEnd) + 1
                        qS = int(start) + 1
                        qE = int(end) + 1
                        ## te_subseq = portion of the seq that matches to TE
                        ## trimmed_seq = portion of the seq that does not match to T
                        #end1 = end + 1
                        te_subseq = seq[start:end+1]
                        if strand == '-':
                            te_subseq = reverse_complement(te_subseq)
                            trimmed_seq = seq[end+1:]
                            trimmed_qual= qual[end+1:]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:start:5' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        else:
                            trimmed_seq = seq[0:start]
                            trimmed_qual= qual[0:start]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:end:5' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        if verbose >= 3: print >> sys.stderr, '1: trimmed: %s %s %s' %(rl_name, trimmed_seq, str(end))
                        if len(trimmed_seq) >= len_cutoff_l:
                            print >> ofile_rr, '%s\t%s\t%s' %(rl_name, tName, strand)
                            print >> ofile_te5, '>%s %s..%s matches %s:%s..%s mismatches:%s\n%s' %(header, qS, qE, TE, tS, tE, mismatch, te_subseq)
                    #query read overlaps 3' end of database TE & trimmed seq > cutoff
                    #elif tEnd == tLen - 1 and (int(start) <= 2 or int(end) >= int(length) - 3) and (length - (match + mismatch)) > len_cutoff and (float(mismatch)/(float(match) + float(mismatch))) <= mismatch_allowance:
                    #elif tEnd >= (tLen - 3) and (int(start) <= 2 or int(end) >= int(length) - 3) and (length - (match + mismatch)) > len_cutoff_m and int(mismatch) <= int(mismatch_allowance):
                    elif tEnd >= (tLen - 3) and (int(start) <= 2 or int(end) >= int(length) - 3) and (match + mismatch) >= len_cutoff_m and (length - (match + mismatch)) >= len_cutoff_l and int (mismatch) <= int(mismatch_allowance):
                        tS = int(tStart) + 1
                        tE = int(tEnd) + 1
                        qS = int(start) + 1
                        qE = int(end) + 1
                        te_subseq = seq[start:end+1]
                        if strand == '-':
                            te_subseq = reverse_complement(te_subseq)
                            trimmed_seq = seq[0:start]
                            trimmed_qual= qual[0:start]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:end:3' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        else:
                            trimmed_seq = seq[end+1:]
                            trimmed_qual= qual[end+1:]
                            seq_id   = header
                            seq_desc = ''
                            seq_id   = '%s:start:3' %(seq_id)
                            header = '%s%s' %(seq_id, seq_desc)
                        if verbose >= 3: print >> sys.stderr, '2: trimmed: %s %s %s' %(rl_name, trimmed_seq, str(end))
                        if len(trimmed_seq) >= len_cutoff_l:
                            print >> ofile_rr, '%s\t%s\t%s' %(rl_name, tName, strand)
                            print >> ofile_te3, '>%s %s..%s matches %s:%s..%s mismatches:%s\n%s' %(header, qS, qE, TE, tS, tE, mismatch, te_subseq)
                    #query read overlaps internal of database TE, no need to trim. These reads pairs will be used as supporting reads
                    #in relocate_align, we get the reads in trimmed files and their pairs.
                    #elif start == 0 and end + 1 == length and (float(mismatch)/(float(match) + float(mismatch))) <= mismatch_allowance:
                    #elif start == 0 and end + 1 == length and int(mismatch) <= int(mismatch_allowance):
                    #for internal match we allowed more mismatch and deletions as the repeat are tend to have mutations among strains
                    elif start <= 2 and end + 1 >= length - 3:
                        trimmed_seq = seq
                        trimmed_qual= qual
                        seq_id   = header
                        seq_desc = ''
                        seq_id   = '%s:middle' %(seq_id)
                        header = '%s%s' %(seq_id, seq_desc)
                        print >> ofile_rr, '%s\t%s\t%s' %(rl_name, tName, strand) 
                        if verbose >= 3: print >> sys.stderr, '3: trimmed: %s %s %s' %(rl_name, trimmed_seq, str(end))
                    ##trimmed reads
                    if len(trimmed_seq) >= len_cutoff_l:
                        print '@%s\n%s\n%s\n%s' %(header, trimmed_seq, qualh, trimmed_qual)
                    ##any read that was in the blat file is written here
                    print >> ofile_fq, '@%s\n%s\n%s\n%s' %(header, seq, qualh, qual)
        filehd.close()
        ofile_te5.close()
        ofile_te3.close()
        ofile_fq.close()
        ofile_rr.close()
if __name__ == '__main__':
    main()

