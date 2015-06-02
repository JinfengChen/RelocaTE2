#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob

def usage():
    test="name"
    message='''
python clean_pairs_memory.py --1 fastq1 --2 fastq2 --repeat path/te_containing_fq --fq_dir original_fq_dir --seqtk pathtoseqtk > unPaired.fq

Takes pairs of fastq files, which is trimmed from repeat blat results seperately, and find their mates in each other, in TE_containing fastqs and original fastqs.
*.matched: contain reads that have mates in eigher trimmed fastq or in original fastq (whether includes these in TE_containing need to test).
*.unPaired.fq: contains trimmed reads that do not have mates found and contain the mate pairs of reads that matched to middle of repeat only if the mate pair is not repeat, but not reads themselve as they are part of repeat (whether includes these mates in TE_containing need to test).
    '''
    print message

#we only need these reads matched to TE but not used as flanking_fq
#store only these reads not labeled with @read_500_403/1:middle or @read_500_403/1:end:5
#in some case read id in pairs are the same, not labeled with /1 or /2, or .r/.f
def parse_fastq(fq_file):
    data =defaultdict(lambda : str)
    #s = re.compile(r':(middle|start|end)')
    s1= re.compile(r'(\S+)(\.[rf])')
    s2= re.compile(r'(\S+)(\/[12])')
    with open (fq_file, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            header = re.split(r' ', line[1:])[0]
            read_t = '1'
            #if s.search(header):
            #    continue
            m1 = s1.search(header)
            m2 = s2.search(header)
            if m1:
                header = m1.groups(0)[0]
                read_t = m1.groups(0)[1]
            elif m2:
                header = m2.groups(0)[0]
                read_t = m2.groups(0)[1]
            #header = re.sub(r'\/.*', '', header)
            seq    = filehd.next().rstrip()
            qualh  = filehd.next().rstrip()
            qual   = filehd.next().rstrip()
            data[header] = read_t
            #print header, read_t
    return data

def parse_fastq_flanking(fq_file):
    s = re.compile(r'(.*):(middle|start|end)')
    s1= re.compile(r'(\S+)\.[rf]')
    s2= re.compile(r'(\S+)\/[12]')
    fq_dict = defaultdict(lambda : list)
    with open (fq_file, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            #header = line[1:]
            #header_to_store = header
            header = re.split(r' ', line[1:])[0]
            m = s.search(header)
            header_to_store = m.groups(0)[0] if m else header
            pos = m.groups(0)[1] if m else 'unknown'
            #print '%s\t%s' %(header_to_store, pos)
            m1 = s1.search(header)
            m2 = s2.search(header)
            if m1:
                header_to_store = m1.groups(0)[0]
            elif m2:
                header_to_store = m2.groups(0)[0]
            #print header_to_store, pos
            #header = re.sub(r'\/.*', '', header)
            seq    = filehd.next().rstrip()
            qualh  = filehd.next().rstrip()
            qual   = filehd.next().rstrip()
            record = '@%s\n%s\n%s\n%s' %(header, seq, qualh, qual)
            fq_dict[header_to_store] = [pos, record, header]
            #print fq_dict[header_to_store][0], fq_dict[header_to_store][1]
    return fq_dict

def parse_fastq_default(fq_file):
    s1= re.compile(r'(\S+)\.[rf]')
    s2= re.compile(r'(\S+)\/[12]')
    fq_dict = defaultdict(lambda : str)
    with open (fq_file, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            #header = line[1:]
            header = re.split(r' ', line[1:])[0]
            header_to_store = header
            m1 = s1.search(header)
            m2 = s2.search(header)
            if m1:
                header_to_store = m1.groups(0)[0]
            elif m2:
                header_to_store = m2.groups(0)[0]
            #header = re.sub(r'\/.*', '', header)
            seq    = filehd.next().rstrip()
            qualh  = filehd.next().rstrip()
            qual   = filehd.next().rstrip()
            record = '@%s\n%s\n%s\n%s' %(header, seq, qualh, qual)
            fq_dict[header_to_store] = record
    return fq_dict


def match_trimmed(fq1, fq2, fq1_match, fq2_match, fq_unPaired, fq_unPaired_info):
    fq1_dict = parse_fastq_flanking(fq1)
    fq2_dict = parse_fastq_flanking(fq2)
    ofile1 = open(fq1_match, 'w')
    ofile2 = open(fq2_match, 'w')
    ofile3 = open(fq_unPaired, 'w')
    ofile4 = open(fq_unPaired_info, 'w')
    #only deal with :end/start here, keep unpaired middle in dictionary
    for hd in sorted(fq2_dict.keys()):
        #print hd
        if not fq2_dict[hd][0] == 'middle':
        #fq2 is not middle
            #print hd
            if fq1_dict.has_key(hd) and fq1_dict[hd][0] != 'middle':
            #paired and both matched to end of repeat, same end or different end?
            ##these informations should be recorded and used as supporting reads
                #print '1'
                print >> ofile1, fq1_dict[hd][1]
                print >> ofile2, fq2_dict[hd][1]
                del fq1_dict[hd]
                del fq2_dict[hd]
            elif fq1_dict.has_key(hd) and fq1_dict[hd][0] == 'middle':
            #paired but mate in fq1 is middle, write fq2 to unpaired and delete both
                #print '2'
                print >> ofile3, fq2_dict[hd][1]
                print >> ofile4, '%s\t%s' %(fq2_dict[hd][2], 2)
                del fq1_dict[hd]
                del fq2_dict[hd]
            else:
            #not paired in trimmed fastq. do nothing, keep id and find pair in orignal fastq
                #print '3'
                pass
        else:
        #fq2 is middle
            if fq1_dict.has_key(hd) and fq1_dict[hd][0] != 'middle':
            #paired but mate in fq2 is middle, write fq1 to unpaired and delete both
            #these informations should be recorded and used as supporting reads
                print >> ofile3, fq1_dict[hd][1]
                print >> ofile4, '%s\t%s' %(fq1_dict[hd][2], 1)
                del fq1_dict[hd]
                del fq2_dict[hd]
            elif fq1_dict.has_key(hd) and fq1_dict[hd][0] == 'middle':
            #paired but both in middle, useless for insertion delete both
                del fq1_dict[hd]
                del fq2_dict[hd]
                pass
            else:
            #not paired in trimmed fastq, do nothing, keep id and find pair in orignal fastq
                pass
    ofile1.close()
    ofile2.close()
    ofile3.close()
    ofile4.close()
    #return dictionary which includes unpaired middles and unpaired end/start
    return (fq1_dict, fq2_dict)

def match_support(fq1_dict, fq2_0, fq2_te, fq1_match, fq2_match, fq_unPaired, fq_unPaired_info, read_flag, fq1_id_temp, fq2_0_temp, seqtk):
    ofile1 = open(fq1_match, 'a')
    ofile2 = open(fq2_match, 'a')
    ofile3 = open(fq_unPaired, 'a')
    ofile5 = open(fq_unPaired_info, 'a')    

    #deal with mates in te_containing fastq 
    fq2_te_dict = parse_fastq(fq2_te)
    for hd in sorted(fq1_dict.keys()):
        if fq1_dict[hd][0] == 'middle':
            if fq2_te_dict.has_key(hd):
            #reads have their mate matched to repeat but not at start/end or middle
            #delete reads, useless for insertion
                del fq1_dict[hd]
        else:
            if fq2_te_dict.has_key(hd):
            #reads have their mate matched to repeat but not at start/end or middle
            #reads matched to start/end, write to unPaired and delete from fq1
                print >> ofile3, fq1_dict[hd][1]
                print >> ofile5, '%s\t%s' %(fq1_dict[hd][2], read_flag)
                del fq1_dict[hd]
    
    #deal with mates in original fastq
    #find the mate of one reads in large fastq file or bam file is slow, how to speed up this. Maybe the only way to speed up is get clipped reads
    #and their mates, unproperly mapped and their mates into fastq and start from there.
    if fq2_te_dict.values()[0] == '1':
        ofile4 = open(fq1_id_temp, 'w')
        for hd in sorted(fq1_dict.keys()):
            print >> ofile4, hd
        ofile4.close()
    else:
        ofile4 = open(fq1_id_temp, 'w')
        for hd in sorted(fq1_dict.keys()):
            print >> ofile4, '%s%s' %(hd, fq2_te_dict.values()[0])
        ofile4.close()
    #get subset of fastq from original fastq
    #cmd = 'seqtk subseq %s %s > %s' %(fq2_0, fq1_id_temp, fq2_0_temp)
    #print cmd
    os.system('%s subseq %s %s > %s' %(seqtk, fq2_0, fq1_id_temp, fq2_0_temp))
    fq2_0_temp_dict = parse_fastq_default(fq2_0_temp)
    #os.system('rm %s %s' %(fq1_id_temp, fq2_0_temp))

    #write paired and unPaired reads 
    for hd in sorted(fq1_dict.keys()):
        if fq1_dict[hd][0] == 'middle':
            if fq2_0_temp_dict.has_key(hd):
            ##use mates as supporting reads
                print >> ofile3, fq2_0_temp_dict[hd]
            else:
            ##no mates found, useless. impossible if input is paired
                pass
        else:
            if fq2_0_temp_dict.has_key(hd):
            ##paired, write both
                print >> ofile1, fq1_dict[hd][1]
                print >> ofile2, fq2_0_temp_dict[hd]
            else:
            ##no mates found, write unpaired. impossible if input is paired
                print >> ofile3, fq1_dict[hd][1]
    ofile1.close()
    ofile2.close()
    ofile3.close()
    ofile5.close()
 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--fq1')
    parser.add_argument('-2', '--fq2')
    parser.add_argument('-r', '--repeat')
    parser.add_argument('-f', '--fq_dir')
    parser.add_argument('-s', '--seqtk')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.fq1) > 0 and len(args.fq2) > 0 and len(args.repeat) and len(args.fq_dir)
    except:
        usage()
        sys.exit(2)

    fastqs = glob.glob('%s/*.f*q*' %(args.fq_dir)) 
    suffix = ''
    s = re.compile(r'(\.f\w*?q.*?$)')
    if s.search(fastqs[0]):
        suffix = s.search(fastqs[0]).groups(0)[0]

    fq1_te = '%s/%s' %(args.repeat, re.sub(r'.flankingReads.fq', r'.ContainingReads.fq', os.path.split(args.fq1)[1]))
    fq2_te = '%s/%s' %(args.repeat, re.sub(r'.flankingReads.fq', r'.ContainingReads.fq', os.path.split(args.fq2)[1]))
    fq1_0  = '%s/%s' %(args.fq_dir, re.sub(r'.te_repeat.flankingReads.fq', r'%s' %(suffix), os.path.split(args.fq1)[1]))
    fq2_0  = '%s/%s' %(args.fq_dir, re.sub(r'.te_repeat.flankingReads.fq', r'%s' %(suffix), os.path.split(args.fq2)[1]))
    fq1_match = '%s.matched' %(args.fq1)
    fq2_match = '%s.matched' %(args.fq2)
    fq_unPaired = '%s.unPaired.fq' %(os.path.splitext(args.fq1)[0])
    fq_unPaired_info = '%s.unPaired.info' %(os.path.splitext(args.fq1)[0])
    #print '%s\n%s\n%s\n%s' %(args.fq2, fq2_te, fq2_0, fq_unPaired)
    
    
    #write pairs that exists in trimmed files (*.flankingReads.fq) to *.matched
    #write pairs that have their mates in trimmed files, but are middles, to *.unPaired.fq
    #Return the id of left reads to get their mates from original fastq
    fq1_dict, fq2_dict = match_trimmed(args.fq1, args.fq2, fq1_match, fq2_match, fq_unPaired, fq_unPaired_info)
    
    #write pairs that have their mates in te_containing fastq in unPaired
    #write pairs get from trimmed and original fastq to *.matched
    fq1_id_temp = '%s.fq1_id_temp.list' %(os.path.splitext(args.fq1)[0])
    fq2_0_temp  = '%s.fq2_0_temp.fq' %(os.path.splitext(args.fq2)[0])
    #print '%s\n%s' %(fq1_id_temp, fq2_0_temp)
    match_support(fq1_dict, fq2_0, fq2_te, fq1_match, fq2_match, fq_unPaired, fq_unPaired_info, 1, fq1_id_temp, fq2_0_temp, args.seqtk)
    fq2_id_temp = '%s.fq2_id_temp.list' %(os.path.splitext(args.fq2)[0])
    fq1_0_temp  = '%s.fq1_0_temp.fq' %(os.path.splitext(args.fq1)[0])
    #print '%s\n%s' %(fq2_id_temp, fq1_0_temp)
    match_support(fq2_dict, fq1_0, fq1_te, fq2_match, fq1_match, fq_unPaired, fq_unPaired_info, 2, fq2_id_temp, fq1_0_temp, args.seqtk)

 
if __name__ == '__main__':
    main()

