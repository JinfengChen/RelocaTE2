#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from collections import OrderedDict
import re
import os
import argparse
import pysam
import glob

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

#Retro1  ACGTC   not.give        Chr4    14199..14203    -       T:7     R:4     L:3     ST:21   SR:9    SL:12
def txt2gff(infile, outfile, ins_type):
    #print infile, outfile
    ofile = open(outfile, 'a')
    count = 0
    r_pos = re.compile(r'(\d+)\.\.(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            #print line
            if len(line) > 2:
                unit = re.split(r'\t',line)
                count += 1
                chro, start, end = ['', 0, 0]
                chro = unit[3]
                strand = unit[5]
                m = r_pos.search(unit[4])
                if m:
                    start = m.groups(0)[0]
                    end   = m.groups(0)[1]
                r_count = re.sub(r'\D+', '', unit[7])
                l_count = re.sub(r'\D+', '', unit[8])
                r_supp  = re.sub(r'\D+', '', unit[10])
                l_supp  = re.sub(r'\D+', '', unit[11])
                r_id    = 'repeat_%s_%s_%s' %(chro, start, end)
                print >> ofile, '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=%s;TSD=%s;Note=%s;Right_junction_reads:%s;Left_junction_reads:%s;Right_support_reads:%s;Left_support_reads:%s;' %(chro, 'RelocaTE2', unit[2], start, end,strand, r_id, unit[1], ins_type, r_count, l_count, r_supp, l_supp)
    ofile.close()

#[name, seq, start, strand]
def Supporting_count(event, tsd_start, teSupportingReads):
    total = 0
    right = 0
    left  = 0
    read1 = []
    read2 = []
    #print 'event: %s' %(event)
    if teSupportingReads.has_key(event):
        #print 'event: %s' %(event)
        for read in teSupportingReads[event]:
            name   = read[0]
            seq    = read[1]
            start  = read[2]
            strand = read[3]
            #print name, seq, start, strand
            if int(start) + len(seq) <= int(tsd_start) and strand == '+':
                total += 1
                left  += 1
                read1.append(name)
            elif int(start) >= int(tsd_start) and strand == '-':
                total += 1
                right += 1
                read2.append(name)
        #del teSupportingReads[event]
        return (total, left, right, ','.join(read1), ','.join(read2))
    else:
        return (0,0,0,'','') 

def read_repeat_name(infiles):
    data = defaultdict(list)
    for infile in infiles:
        with open (infile, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2: 
                    unit = re.split(r'\t',line)
                    data[unit[0]] = [unit[1], unit[2]]
    return data

#s1= re.compile(r'(\S+)\.[rf]')
#s2= re.compile(r'(\S+)\/[12]')
## define insertion as repeat family according reads <-> repeats relation from blat or bam files
def insertion_family_supporting(reads, read_repeat):
    repeat_family = defaultdict(lambda : int())
    for read in re.split(r',', reads):
        read_name  = read
        read_name1 = '%s/1' %(read)
        read_name2 = '%s/2' %(read)
        read_name3 = '%s.f' %(read)
        read_name4 = '%s.r' %(read)
        if read_repeat.has_key(read_name):
            repeat_family[read_repeat[read_name][0]] += 1
        elif read_repeat.has_key(read_name1):
            repeat_family[read_repeat[read_name1][0]] += 1
            #print '%s,%s,%s,' %(read_name, read_name1, read_repeat[read_name1])
        elif read_repeat.has_key(read_name2):
            repeat_family[read_repeat[read_name2][0]] += 1
            #print '%s,%s,%s,' %(read_name, read_name2, read_repeat[read_name2])
        elif read_repeat.has_key(read_name3):
            repeat_family[read_repeat[read_name3][0]] += 1
        elif read_repeat.has_key(read_name4):
            repeat_family[read_repeat[read_name4][0]] += 1
    if len(repeat_family.keys()) == 1:
        #return first element if only have one repeat
        return repeat_family.keys()[0]
    elif len(repeat_family.keys()) > 1:
        #return the one have largest value
        sorted_by_value = OrderedDict(sorted(repeat_family.items(), key=lambda x: x[1]))
        return sorted_by_value.keys()[-1]
    else:
        return 'NA'


## define insertion as repeat family according reads <-> repeats relation from blat or bam files
def insertion_family(reads, read_repeat):
    repeat_family = defaultdict(lambda : int())
    r = re.compile(r'(.*):(start|end):(5|3)')
    for read in re.split(r',', reads):
        m = r.search(read)
        read_name = m.groups(0)[0] if m else 'NA'
        if read_name != 'NA' and read_repeat.has_key(read_name):
            repeat_family[read_repeat[read_name][0]] += 1
    if len(repeat_family.keys()) == 1:
        #return first element if only have one repeat
        return repeat_family.keys()[0]
    elif len(repeat_family.keys()) > 1:
        #return the one have largest value
        sorted_by_value = OrderedDict(sorted(repeat_family.items(), key=lambda x: x[1]))
        return sorted_by_value.keys()[-1]
    else:
        return ''

#mPing   GAA     not.give        Chr4    3386246..3386248        -       T:3     R:1     L:2     ST:0    SR:0    SL:0
def write_output(top_dir, result, read_repeat, usr_target, exper, TE, required_reads, required_left_reads, required_right_reads, teInsertions, teInsertions_reads, teSupportingReads, existingTE_inf, existingTE_found, teReadClusters, bedtools, lib_size):
    ref     = '%s/%s.%s.all_ref_insert.txt' %(result, usr_target, TE)
    ref_gff = '%s/%s.%s.all_ref_insert.gff' %(result, usr_target, TE)
    REF     = open ('%s/%s.%s.all_ref_insert.txt' %(result, usr_target, TE), 'w')
    #REFGFF     = open ('%s/%s.%s.all_ref_insert.gff' %(result, usr_target, TE), 'w')
    #READS      = open ('%s/%s.%s.all_ref_reads.list' %(result, usr_target, TE), 'w')
    r = re.compile(r'(\w+):(\d+)-(\d+):(.*)')
    for te_id in sorted(existingTE_found.keys()):
        #print 'Found: %s' %(te_id)
        ##junction reads
        strand, chro, start, end = ['', '', '', '']
        if r.search(te_id):
            strand = r.search(te_id).groups(0)[3]
            chro   = r.search(te_id).groups(0)[0]
            start  = r.search(te_id).groups(0)[1]
            end    = r.search(te_id).groups(0)[2]
        #print '%s\t%s\t%s\t%s' %(strand, chro, start, end)
        l_count, r_count = [0, 0]
        if strand == '+':
            l_count = len(existingTE_found[te_id]['start'].keys())
            r_count = len(existingTE_found[te_id]['end'].keys())
        else:
            l_count = len(existingTE_found[te_id]['end'].keys())
            r_count = len(existingTE_found[te_id]['start'].keys())
        #print '%s\t%s' %(l_count, r_count)
        ##reads and events
        total_supporting_l, left_supporting_l, right_supporting_l, left_reads_l, right_reads_l = [0, 0, 0, '', '']
        total_supporting_r, left_supporting_r, right_supporting_r, left_reads_r, right_reads_r = [0, 0, 0, '', '']
        reads   = defaultdict(lambda : int())
        event_l = defaultdict(lambda : int())
        event_r = defaultdict(lambda : int())
        if len(existingTE_found[te_id]['start'].keys()) > 0:
            for rd in existingTE_found[te_id]['start'].keys():
                reads[rd] = 0
                event_l[existingTE_found[te_id]['start'][rd]] += 1
            event_l_sorted = OrderedDict(sorted(event_l.items(), key=lambda x:x[1]))
            event_l_top    = event_l_sorted.keys()[-1]
            total_supporting_l, left_supporting_l, right_supporting_l, left_reads_l, right_reads_l = Supporting_count(event_l_top, start, teSupportingReads)
        if len(existingTE_found[te_id]['end'].keys()) > 0:
            for rd in existingTE_found[te_id]['end'].keys():
                reads[rd] = 1
                event_r[existingTE_found[te_id]['end'][rd]] += 1
            event_r_sorted = OrderedDict(sorted(event_r.items(), key=lambda x:x[1]))
            event_r_top    = event_r_sorted.keys()[-1]
            total_supporting_r, left_supporting_r, right_supporting_r, left_reads_r, right_reads_r = Supporting_count(event_r_top, str(int(end)+1), teSupportingReads)
        #event_l_sorted = OrderedDict(sorted(event_l.items(), key=lambda x:x[1]))
        #event_r_sorted = OrderedDict(sorted(event_r.items(), key=lambda x:x[1]))
        #event_l_top    = event_l_sorted.keys()[-1]
        #event_r_top    = event_r_sorted.keys()[-1]
        ##repeat family
        repeat_junction = insertion_family(','.join(reads.keys()), read_repeat)
        #print 'repeat_junction: %s' %(repeat_junction)
        ##supporting reads
        #total_supporting_l, left_supporting_l, right_supporting_l, left_reads_l, right_reads_l = Supporting_count(event_l_top, start, teSupportingReads)
        #total_supporting_r, left_supporting_r, right_supporting_r, left_reads_r, right_reads_r = Supporting_count(event_r_top, str(int(end)+1), teSupportingReads)
        #print 'left: t:%s\tl:%s\tr:%s\tl:%s\tr:%s' %(total_supporting_l, left_supporting_l, right_supporting_l, left_reads_l, right_reads_l)
        #print 'right: t:%s\tl:%s\tr:%s\tl:%s\tr:%s' %(total_supporting_r, left_supporting_r, right_supporting_r, left_reads_r, right_reads_r)
        ##output
        if l_count > 0 and r_count > 0:
            print >> REF, '%s\t%s\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_junction, 'TSD', exper, chro, start, end, strand, l_count+r_count, r_count, l_count, left_supporting_l+right_supporting_r, right_supporting_r, left_supporting_l)
    REF.close()
    txt2gff(ref, ref_gff, 'Shared, in ref and reads') 
 
def read_direction(strands):
    plus  = 0
    minus = 0
    for s in strands:
        if s == '+':
            plus += 1
        else:
            minus += 1
    if plus > 0 and plus > minus:
        return 'left'
    elif minus > 0 and minus > plus:
        return 'right'

def get_boundary(reads_list, direction):
    read_starts = list()
    read_ends   = list()
    for read_inf in reads_list:
        read_starts.append(int(read_inf[2]))
        read_ends.append(int(read_inf[2])+len(read_inf[1]))
        #print '%s\t%s\t%s\t%s' %(read_inf[0], read_inf[1], read_inf[2], read_inf[3])
    if direction == 'right':
        #print 'right: %s' %(sorted(read_starts, key=int)[0])
        return sorted(read_starts, key=int)[0]
    elif direction == 'left':
        #print 'left: %s' %(sorted(read_ends, key=int)[-1])
        return sorted(read_ends, key=int)[-1]

def read_to_repeat(read_name, read_repeat):
    read_name1 = '%s/1' %(read_name)
    read_name2 = '%s/2' %(read_name)
    read_name3 = '%s.f' %(read_name)
    read_name4 = '%s.r' %(read_name)
    repeat_family = []
    if read_repeat.has_key(read_name):
        repeat_family = read_repeat[read_name]
    elif read_repeat.has_key(read_name1):
        repeat_family = read_repeat[read_name1]
    elif read_repeat.has_key(read_name2):
        repeat_family = read_repeat[read_name2]
    elif read_repeat.has_key(read_name3):
        repeat_family = read_repeat[read_name3]
    elif read_repeat.has_key(read_name4):
        repeat_family = read_repeat[read_name4]
    return repeat_family

def TSD_from_read_depth(r, read_repeat, teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found):
    #determine TSD from read depth at insertions site
    #count depth to find TSD in
    #if there are 5 reads (2 right, 3 left) they
    #should only be a depth of 5 at the TSD
    #teReadCluster = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str()))))
    #teReadClusters[event]['read_inf'][name]['strand']= strand
    #teReadClusters_count[event]['read_count'] += 1
    #teReadClusters_depth[event]['read_inf']['depth'][i] += 1

    #check how many insertion in each cluster by find start/end of junction reads
    #split cluster into subcluster, then find tsd using depth method
    #reestimate supporting reads for each cluster
    r5 = re.compile(r'start:[53]$')
    r3 = re.compile(r'end:[53]$')
 
    for cluster in sorted(teReadClusters.keys(), key=int):
        chro    = teReadClusters[cluster]['read_inf']['seq']['chr']
        ##check insertion number
        left_reads = defaultdict(lambda : list())
        right_reads= defaultdict(lambda : list())
        
        teReadClusters_sub = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str()))))
        teReadClusters_sub_count = defaultdict(lambda : defaultdict(lambda : int()))
        teReadClusters_sub_depth = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int()))))
        teReadClusters_sub_type  = defaultdict(lambda : int())
        #print 'tsdfiner%s' %(cluster)
        for name in teReadClusters[cluster]['read_inf'].keys():
            if name == 'seq':
                #skip empty line
                continue
            seq    = teReadClusters[cluster]['read_inf'][name]['seq']
            start  = teReadClusters[cluster]['read_inf'][name]['start']
            strand = teReadClusters[cluster]['read_inf'][name]['strand']
            #print name, start, seq, strand
            end    = int(start) + len(seq)
            if strand == '+':
                if r5.search(name):
                    pos  = 'right'
                    right_reads[start].append(name)
                elif r3.search(name):
                    pos  = 'left'
                    left_reads[end].append(name)
            elif strand == '-':
                if r5.search(name):
                    pos  = 'left'
                    left_reads[end].append(name)
                elif r3.search(name):
                    pos  = 'right'
                    right_reads[start].append(name)
        #print len(left_reads.keys()), len(right_reads.keys())
        if (len(left_reads.keys()) > 1 and len(right_reads.keys()) >= 1) or (len(left_reads.keys()) >= 1 and len(right_reads.keys()) > 1):
            ##more than two junction
            count_tsd = 0
            pairs_tsd = defaultdict(lambda : int())
            ##find pairs for one insertions
            for start1 in left_reads.keys():
                min_dist = 0
                min_pair = ''
                for start2 in right_reads.keys():
                    if min_dist == 0: 
                        min_dist = abs(int(start2) - int(start1))
                        min_pair = start2
                    elif min_dist > abs(int(start2) - int(start1)):
                        min_dist = abs(int(start2) - int(start1))
                        min_pair = start2
                if min_dist <= 100:
                    ##find pairs
                    count_tsd += 1
                    pairs_tsd[start1]   = 1
                    pairs_tsd[min_pair] = 1
                    teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 2
                    for read in left_reads[start1]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                    for read in right_reads[min_pair]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                else:
                    ##do not find pairs
                    count_tsd += 1
                    pairs_tsd[start1]   = 1
                    teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                    for read in left_reads[start1]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
            ##set unpaired
            for start2 in right_reads.keys():
                if not pairs_tsd.has_key(start2):
                    #not paired junction
                    count_tsd += 1
                    teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                    for read in right_reads[start2]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
        elif len(left_reads.keys()) > 1:
            count_tsd = 0
            for start1 in left_reads.keys():
                count_tsd += 1
                teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
        elif len(right_reads.keys()) > 1:
            count_tsd = 0
            for start2 in right_reads.keys():
                count_tsd += 1
                teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
        elif len(left_reads.keys()) == 1 and len(right_reads.keys()) == 1:
            ##one right and one left junction
            #print left_reads.keys()[0], right_reads.keys()[0]
            if abs(int(left_reads.keys()[0]) - int(right_reads.keys()[0])) > 100:
                ##two far from each other, might be one end from two insertion
                teReadClusters_sub_type['%s-1' %(cluster)] = 1
                teReadClusters_sub_type['%s-2' %(cluster)] = 1
                start1 = left_reads.keys()[0]
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                start2 = right_reads.keys()[0]
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-2' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-2' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-2' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-2' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
            else:
                ##one junction
                teReadClusters_sub_type['%s-1' %(cluster)] = 2
                start1 = left_reads.keys()[0]
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                start2 = right_reads.keys()[0]
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth) 
        else:
            ##one junction with one end support
            if len(left_reads.keys()) > 0:
                teReadClusters_sub_type['%s-1' %(cluster)] = 1
                start1 = left_reads.keys()[0] 
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
            elif len(right_reads.keys()) > 0:
                start2 = right_reads.keys()[0]
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)

        #print 'TSD finder: %s' %(cluster)
        ###teReadCluster_sub_depth add above
        ###deal with teReadClusters_sub, still store at cluster in teInsertions, write a TSD_check for this only.
        for sub_cluster in teReadClusters_sub_depth.keys():
            #print sub_cluster
            
            #TSD_len = 0
            #read_total = teReadClusters_sub_count[sub_cluster]['read_count']
            #for chrs_pos in sorted(teReadClusters_sub_depth[sub_cluster]['read_inf']['depth'].keys(), key=int):
            #    depth = teReadClusters_sub_depth[sub_cluster]['read_inf']['depth'][chrs_pos]
            #    if float(depth) >= 0.6*float(read_total):
            #        TSD_len += 1
            if teReadClusters_sub_type[sub_cluster] == 1:
                TSD = 'UKN'
                for name in teReadClusters_sub[sub_cluster]['read_inf'].keys():
                    real_name = r.search(name).groups(0)[0] if r.search(name) else ''
                    seq    = teReadClusters_sub[sub_cluster]['read_inf'][name]['seq']
                    start  = teReadClusters_sub[sub_cluster]['read_inf'][name]['start']
                    strand = teReadClusters_sub[sub_cluster]['read_inf'][name]['strand']
                    TSD_check_single(cluster, seq, chro, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            else:
                TSD_len = 0
                if tsd_finder(sub_cluster, 1, teReadClusters_sub_count, teReadClusters_sub_depth):
                    TSD_len = tsd_finder(sub_cluster, 1, teReadClusters_sub_count, teReadClusters_sub_depth)
                elif tsd_finder(sub_cluster, 0.8, teReadClusters_sub_count, teReadClusters_sub_depth):
                    TSD_len = tsd_finder(sub_cluster, 0.8, teReadClusters_sub_count, teReadClusters_sub_depth)
                elif tsd_finder(sub_cluster, 0.6, teReadClusters_sub_count, teReadClusters_sub_depth):
                    TSD_len = tsd_finder(sub_cluster, 0.6, teReadClusters_sub_count, teReadClusters_sub_depth)

                if TSD_len > 0:
                    #print TSD_len
                    TSD = '.'*TSD_len
                    for name1 in teReadClusters_sub[sub_cluster]['read_inf'].keys():
                        real_name = r.search(name1).groups(0)[0] if r.search(name1) else ''
                        seq    = teReadClusters_sub[sub_cluster]['read_inf'][name1]['seq']
                        start  = teReadClusters_sub[sub_cluster]['read_inf'][name1]['start']
                        strand = teReadClusters_sub[sub_cluster]['read_inf'][name1]['strand']
                        #print name1, seq, start, strand, chro
                        TSD_check(cluster, seq, chro, start, real_name, read_repeat, name1, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
                        #print 'after tsd_check'
                else:          
                    #what if we can not find TSD? still could be insertions
                    TSD = 'UKN'
                    #print TSD
                    for name in teReadClusters_sub[sub_cluster]['read_inf'].keys():
                        real_name = r.search(name).groups(0)[0] if r.search(name) else ''
                        seq    = teReadClusters_sub[sub_cluster]['read_inf'][name]['seq']
                        start  = teReadClusters_sub[sub_cluster]['read_inf'][name]['start']
                        strand = teReadClusters_sub[sub_cluster]['read_inf'][name]['strand']
                        TSD_check(cluster, seq, chro, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            #print 'End of cycle'

def tsd_finder(sub_cluster, tsd_depth, teReadClusters_sub_count, teReadClusters_sub_depth):
    TSD_len  = 0
    read_total = teReadClusters_sub_count[sub_cluster]['read_count']
    for chrs_pos in sorted(teReadClusters_sub_depth[sub_cluster]['read_inf']['depth'].keys(), key=int):
        depth = teReadClusters_sub_depth[sub_cluster]['read_inf']['depth'][chrs_pos]
        if float(depth) >= float(tsd_depth)*float(read_total):
            TSD_len += 1
    return TSD_len

def align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, chro, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads):
    range_allowance = 1000
    padded_start    = bin_ins[0] - range_allowance
    padded_end      = bin_ins[-1] + range_allowance 
    #insertions
    #print 'insertions: %s' %(name)
    if (int(start) >= padded_start and int(start) <= padded_end) or (int(end) >= padded_start and int(end) <= padded_end):
        bin_ins.extend([int(start), int(end)])
        bin_ins = sorted(bin_ins, key=int)
        if r.search(name):
            real_name = r.search(name).groups(0)[0]
            if not r_tsd.search(TSD):
                TSD_check(count, seq, chro, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            else:
                calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
        elif not r.search(name) and not record.is_paired:
            #reads not matched to repeat and not mates of junctions
            #reads are mates of reads matched to middle of repeat
            #supporting reads
            teSupportingReads[count].append([name, seq, start, strand])
    else:
        #if start and end do not fall within last start and end
        #we now have a different insertion event
        count += 1
        if r.search(name):
            real_name = r.search(name).groups(0)[0]
            if not r_tsd.search(TSD):
                TSD_check(count, seq, chro, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            else:
                calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
        elif not r.search(name) and not record.is_paired:
            #reads not matched to repeat and not mates of junctions
            #reads are mates of reads matched to middle of repeat
            #supporting reads
            teSupportingReads[count].append([name, seq, start, strand])
        #initial insertion site boundary
        bin_ins = [int(start), int(end)]
        #print '%s\t%s' %(count, bin_ins)
    return (bin_ins, count)

#existing_TE_bed_reader
#Chr4    1072    1479    Simple_repeat:1072-1479 1       +
#Chr4    1573    1779    Simple_repeat:1573-1779 0       +
def existing_TE_bed_reader(infile, existingTE_intact, chro):
    #print 'Reading existing TE bed'
    with open(infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip() 
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                ##intact repeat
                #print line
                if int(unit[4]) == 1 and unit[0] == chro:
                    te_id = '%s:%s-%s:%s' %(unit[0], unit[1], unit[2], unit[5])
                    existingTE_intact[unit[0]]['start'][int(unit[1])] = te_id
                    existingTE_intact[unit[0]]['end'][int(unit[2])] = te_id
                    #print '%s\t%s\t%s' %(unit[0], unit[1], 'start')
                    #print '%s\t%s\t%s' %(unit[0], unit[2], 'end')
                    #if unit[5] == '+':
                    #    existingTE_intact[unit[0]]['start'][unit[1]] = te_id
                    #    existingTE_intact[unit[0]]['end'][unit[2]] = te_id
                    #else:
                    #    existingTE_intact[unit[0]]['end'][unit[1]] = te_id
                    #    existingTE_intact[unit[0]]['start'][unit[2]] = te_id

def existingTE_RM_ALL(top_dir, infile, existingTE_inf):
    ofile_RM = open('%s/existingTE.bed' %(top_dir), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\s+',line)
                #print line
                #print unit[5], unit[9], unit[12], unit[13], unit[14]
                if unit[9] == '+':
                    for i in range(int(unit[6])-2, int(unit[6])+3):
                        existingTE_inf[unit[5]]['start'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7], '1', '+')
                    #print unit[10], 'start', unit[6]
                    for i in range(int(unit[7])-2, int(unit[7])+3):
                        existingTE_inf[unit[5]]['end'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '+')
                    #print unit[10], 'end', unit[7]
                    print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])), str(int(unit[7])), unit[11],unit[6],unit[7], '1', '+')
                elif unit[9] == 'C':
                    for i in range(int(unit[6])-2, int(unit[6])+3):
                        existingTE_inf[unit[5]]['start'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7],'1', '-')
                    #print unit[10], 'start', unit[6]
                    for i in range(int(unit[7])-2, int(unit[7])+3):
                        existingTE_inf[unit[5]]['end'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '-')
                    #print unit[10], 'end', unit[7]
                    print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])), str(int(unit[7])), unit[11],unit[6],unit[7], '1', '-')
    ofile_RM.close()

def existingTE_RM(top_dir, infile, existingTE_inf):
    r_end = re.compile(r'\((\d+)\)')
    ofile_RM = open('%s/existingTE.bed' %(top_dir), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\s+',line)
                #print line
                #print unit[5], unit[9], unit[12], unit[13], unit[14]
                if unit[9] == '+':
                    if int(unit[12]) == 1:
                        for i in range(int(unit[6])-2, int(unit[6])+3):
                            existingTE_inf[unit[5]]['start'][int(i)] = 1
                        print >> ofile_RM, '%s\t%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], 'start', str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7], '1', '+')
                        #print unit[10], 'start', unit[6]
                    if len(unit[14]) == 3:
                        unit[14] =re.sub(r'\(|\)', '', unit[14])
                        if int(unit[14]) == 0:
                            for i in range(int(unit[7])-2, int(unit[7])+3):
                                existingTE_inf[unit[5]]['end'][int(i)] = 1
                            print >> ofile_RM, '%s\t%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], 'end', str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '+')
                            #print unit[10], 'end', unit[7]
                elif unit[9] == 'C':
                    if len(unit[12]) == 3:
                        unit[12] =re.sub(r'\(|\)', '', unit[12])
                        if int(unit[12]) == 0:
                            for i in range(int(unit[6])-2, int(unit[6])+3):
                                existingTE_inf[unit[5]]['start'][int(i)] = 1
                            print >> ofile_RM, '%s\t%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], 'start', str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7],'1', '-')
                            #print unit[10], 'start', unit[6]
                    if int(unit[14]) == 1:
                        for i in range(int(unit[7])-2, int(unit[7])+3):
                            existingTE_inf[unit[5]]['end'][int(i)] = 1
                        print >> ofile_RM, '%s\t%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], 'end', str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '-')
                        #print unit[10], 'end', unit[7]
    ofile_RM.close() 

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    for i in range(len(bases)):
        bases[i] = complement[bases[i]] if complement.has_key(bases[i]) else bases[i]
    return ''.join(bases)

def reverse_complement(seq):
    return complement(seq[::-1])

def calculate_cluster_depth(event, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth):
    teReadClusters_count[event]['read_count'] += 1
    teReadClusters[event]['read_inf'][name]['seq']   = seq
    teReadClusters[event]['read_inf'][name]['start'] = start
    teReadClusters[event]['read_inf'][name]['strand']= strand
    for i in range(int(start), int(start)+len(seq)):
        teReadClusters_depth[event]['read_inf']['depth'][i] += 1

def TSD_check_single(event, seq, chro, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found):
    ##TSD already specified by usr, set tsd for single end of junction. actural the end of read if read on left of junction or start if read on right.
    ##seq is entire trimmd read, not just the TSD portion of the read
    ##start is the first postition of the entire read match to ref
    repeat = read_repeat[real_name] # need to deal with any te, get infor from name of reads
    rev_com = reverse_complement(seq)
    result    = 0
    pos       = ''
    TE_orient = 0
    TSD_start = 0
    TSD_seq   = ''
    r5 = re.compile(r'start:[53]$')
    r3 = re.compile(r'end:[53]$')
    #r5_tsd = re.compile(r'^(%s)' %(TSD))
    #r3_tsd = re.compile(r'(%s)$' %(TSD))
    #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    ##start means that the TE was removed from the start of the read
    ##5 means the trimmed end mapps to the 5prime end of the TE
    ##3 means the trimmed end mapps to the 3prime end of the TE
    if strand == '+':
        if r5.search(name):
            result    = 1
            TSD_seq   = 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
        elif r3.search(name):
            result    = 1
            TSD_seq   = 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            TSD_start = int(start) + len(seq)
    elif strand == '-':
        if r5.search(name):
            result    = 1
            TSD_seq   = 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            TSD_start = int(start) + len(seq)
        elif r3.search(name):
            result    = 1
            TSD_seq   = 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
    #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    if result and TE_orient:
        tir1_end, tir2_end = [0, 0]
        #if 0:
        #    continue
        if pos == 'left':
            tir1_end = int(TSD_start)
            #print 'tir1: %s' %(tir1_end)
        elif pos == 'right':
            tir2_end = int(TSD_start) - 1
            #print 'tir2: %s' %(tir2_end)
        if tir1_end > 0 and existingTE_inf[chro]['start'].has_key(tir1_end):
            te_id = existingTE_inf[chro]['start'][tir1_end]
            existingTE_found[te_id]['start'][name] = event
            #print 'tir1'
        elif tir2_end > 0 and existingTE_inf[chro]['end'].has_key(tir2_end):
            te_id = existingTE_inf[chro]['end'][tir2_end]
            existingTE_found[te_id]['end'][name] = event
            #print 'tir2'
        else:
            #print 'not match'
            ##non reference insertions
            teInsertions[event][TSD_start][TSD_seq]['count']   += 1   ## total junction reads
            teInsertions[event][TSD_start][TSD_seq][pos]       += 1   ## right/left junction reads
            teInsertions[event][TSD_start][TSD_seq][TE_orient] += 1   ## plus/reverse insertions
            #read_name = re.sub(r':start|:end', '', name)
            teInsertions_reads[event][TSD_start][TSD_seq]['read'].append(name)
            #print '1: %s\t 2: %s' %(read_name, teInsertions_reads[event][TSD_seq][TSD_start]['read'])
            #print 'C: %s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient)

def TSD_check(event, seq, chro, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found):
    ##TSD already specified by usr, not unknown
    ##seq is entire trimmd read, not just the TSD portion of the read
    ##start is the first postition of the entire read match to ref
    repeat = read_repeat[real_name] # need to deal with any te, get infor from name of reads
    rev_com = reverse_complement(seq)
    result    = 0
    pos       = ''
    TE_orient = 0
    TSD_start = 0
    TSD_seq   = ''
    r5 = re.compile(r'start:[53]$')
    r3 = re.compile(r'end:[53]$')
    r5_tsd = re.compile(r'^(%s)' %(TSD))
    r3_tsd = re.compile(r'(%s)$' %(TSD))
    #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    ##start means that the TE was removed from the start of the read
    ##5 means the trimmed end mapps to the 5prime end of the TE
    ##3 means the trimmed end mapps to the 3prime end of the TE
    if strand == '+':
        if r5.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
            result    = 1
            TSD_seq   = r5_tsd.search(seq).groups(0)[0] if r5_tsd.search(seq) else 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
        elif r3.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
            result    = 1
            TSD_seq   = r3_tsd.search(seq).groups(0)[0] if r3_tsd.search(seq) else 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            TSD_start = int(start) + (len(seq)-len(TSD))
    elif strand == '-':
        if r5.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
            result    = 1
            TSD_seq   = r3_tsd.search(seq).groups(0)[0] if r3_tsd.search(seq) else 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            TSD_start = int(start) + (len(seq)-len(TSD))
        elif r3.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
            result    = 1
            TSD_seq   = r5_tsd.search(seq).groups(0)[0] if r5_tsd.search(seq) else 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
    #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    if result and TE_orient:
        tir1_end, tir2_end = [0, 0]
        #if 0:
        #    continue
        if pos == 'left':
            tir1_end = int(start) + len(seq)
            #print 'tir1: %s' %(tir1_end)
        elif pos == 'right':
            tir2_end = int(start) - 1
            #print 'tir2: %s' %(tir2_end)
        if tir1_end > 0 and existingTE_inf[chro]['start'].has_key(tir1_end):
            te_id = existingTE_inf[chro]['start'][tir1_end]
            existingTE_found[te_id]['start'][name] = event
            #print 'tir1'
        elif tir2_end > 0 and existingTE_inf[chro]['end'].has_key(tir2_end):
            te_id = existingTE_inf[chro]['end'][tir2_end]
            existingTE_found[te_id]['end'][name] = event 
            #print 'tir2'
        else:
            #print 'not match'
            ##non reference insertions
            teInsertions[event][TSD_start][TSD_seq]['count']   += 1   ## total junction reads
            teInsertions[event][TSD_start][TSD_seq][pos]       += 1   ## right/left junction reads
            teInsertions[event][TSD_start][TSD_seq][TE_orient] += 1   ## plus/reverse insertions
            #read_name = re.sub(r':start|:end', '', name)
            teInsertions_reads[event][TSD_start][TSD_seq]['read'].append(name)
            #print '1: %s\t 2: %s' %(read_name, teInsertions_reads[event][TSD_seq][TSD_start]['read'])
            #print 'C: %s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient)

def convert_tag(tag):
    tags = {}
    for t in tag:
        tags[t[0]] = t[1]
    return tags

def find_insertion_cluster_bam(align_file, read_repeat, target, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found, teSupportingReads):
    r = re.compile(r'(.*):(start|end):(5|3)')
    r_tsd = re.compile(r'UNK|UKN|unknown', re.IGNORECASE)
    r_cg  = re.compile(r'[SID]')
    bin_ins        = [0]
    count          = 0
    TSD_len        = len(TSD)

    ref  = 'None' if target == 'ALL' else target
    fsam = pysam.AlignmentFile(align_file, 'rb')
    rnames = fsam.references
    for record in fsam.fetch(reference=ref, until_eof = True):
        if not record.is_unmapped:
            name   = record.query_name
            flag   = record.flag
            start  = int(record.reference_start) + 1
            MAPQ   = record.mapping_quality
            cigar  = record.cigarstring
            seq    = record.query_sequence
            tag    = record.tags if record.tags else []
            length = len(seq)
            chro   = rnames[record.reference_id]
            end    = int(start) + int(length) - 1 #should not allowed for indel or softclip
            strand = ''
            # flag is 0 is read if read is unpaired and mapped to plus strand
            if int(flag) == 0:
                strand = '+'
            else:
                strand = '-' if record.is_reverse else '+'
            if r_cg.search(cigar):
                continue
            tags = convert_tag(tag)
            #print '%s\t%s\t%s' %(name, start, length)
            # filter low quality mapping reads: 
            # 1. paired-end reads at least have one reads unique mapped (MAPQ set to 0 for both reads if both are repeat, else should be > 0 at least one unique mapped)
            # 2. unpaired reads should unique mapped, no gap, mismatch <= 3 and no suboptimal alignment
            #print 'before: %s\t%s\t%s' %(name, count, bin_ins)
            #if record.is_proper_pair and (int(MAPQ) >= 29 or tags['XT'] == 'U'):
            if record.is_proper_pair and int(MAPQ) > 0:
                bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, chro, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
            elif not record.is_paired:
                #if tags['XT'] == 'U' and int(tags['XO']) == 0 and int(tags['XM']) <= 3 and int(tags['X1']) == 0:
                #if tags['XT'] == 'U' and int(tags['XO']) == 0 and (int(tags['XM']) <= 3 or int(tags['X1']) == 0):
                if tags['XT'] == 'U' and int(tags['XO']) == 0 and int(tags['X1']) <= 3:
                #if tags['XT'] == 'U':
                    bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, chro, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
            teReadClusters[count]['read_inf']['seq']['chr'] = chro
            #print 'after: %s\t%s\t%s' %(name, count, bin_ins)

    ###TSD not given we infer from read depth
    if r_tsd.search(TSD):
        TSD_from_read_depth(r, read_repeat, teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)        
        
def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

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


def main():

    required_reads       = 1           ## rightReads + leftReads needs to be > to this value
    required_left_reads  = 1           ## needs to be >= to this value
    required_right_reads = 1           ## needs to be >= to this value
    align_file           = sys.argv[1] ## combined bowtie or bwa results, sam format only 
    usr_target           = sys.argv[2] ## chromosome to analyze: ALL or Chr1..N
    genome_path          = sys.argv[3] ## genome sequence
    TE                   = sys.argv[4] ## repeat to analyze: ALL or mPing/other te name 
    regex_file           = sys.argv[5] ## regex.txt
    exper                = sys.argv[6] ## prefix for output, title: HEG4
    flank_len            = sys.argv[7] ## length of seq flanking insertions to be returned: 100
    existing_TE          = sys.argv[8] ## existingTE.blatout
    mm_allow             = sys.argv[9] ## mismatches allowed: 0, 1, 2, 3
    bowtie2              = sys.argv[10] ## use bowtie2 or not: 1 or 0
    lib_size             = sys.argv[11] ## insert size of library
    #relax_reference      = sys.argv[11]## relax mode for existing TE: 1 or 0
    #relax_align          = sys.argv[12]## relax mode for insertion: 1 or 0
    bowtie_sam           = 1           ## change to shift or remove in V2
    existingTE_inf       = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str)))
    existingTE_found     = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str)))
    bwa                  = 0
    teInsertions         = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int()))))
    teInsertions_reads   = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list()))))
    teReadClusters       = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str()))))
    teReadClusters_count = defaultdict(lambda : defaultdict(lambda : int()))
    teReadClusters_depth = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int()))))
    teSupportingReads    = defaultdict(lambda : list())

    top_dir = re.split(r'/', os.path.dirname(os.path.abspath(align_file)))[:-1]
    #read existing TE from file
    #existing_TE_intact = defaultdict(lambda : defaultdict(lambda : int()))
    existingTE_intact = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str)))
    existingTE_bed     = '%s/existingTE.bed' %('/'.join(top_dir))
    existingTE_bed_chr = '%s/existingTE.%s.bed' %('/'.join(top_dir), usr_target)
    #print existingTE_bed
    if os.path.isfile(existingTE_bed_chr) and os.path.getsize(existingTE_bed_chr) > 0:
        existing_TE_bed_reader(existingTE_bed_chr, existingTE_intact, usr_target)
    else:
        os.system('grep -P \"%s\\t\" %s > %s' %(usr_target, existingTE_bed, existingTE_bed_chr))
        existing_TE_bed_reader(existingTE_bed_chr, existingTE_intact, usr_target)
        os.system('rm %s' %(existingTE_bed_chr))
        #print 'Existing TE file does not exists or zero size'

    bedtools = ''
    try:
        subprocess.check_output('which bedtools', shell=True)
        bedtools = subprocess.check_output('which bedtools', shell=True)
        bedtools = re.sub(r'\n', '', bedtools)
    except:
        bedtools = '/opt/bedtools/2.17.0-25-g7b42b3b/bin//bedtools'


    ##get the regelar expression patterns for mates and for the TE
    ##when passed on the command line as an argument, even in single
    ##quotes I lose special regex characters
    s         = re.compile(r'[\[.*+?]')
    mate_file = parse_regex(regex_file)
    TSD       = mate_file[3]
    TSDpattern= 1 if s.search(TSD) else 0


    ##read -> repeat relation
    #top_dir = re.split(r'/', os.path.dirname(os.path.abspath(align_file)))[:-1]
    result  = '%s/results' %('/'.join(top_dir))
    read_repeat_files = []
    if usr_target == 'ALL':
        read_repeat_files = glob.glob('%s/te_containing_fq/*.read_repeat_name.split.txt' %('/'.join(top_dir)))
    else:
        read_repeat_files = glob.glob('%s/te_containing_fq/%s.read_repeat_name.split.txt' %('/'.join(top_dir), usr_target))
    #read_repeat_files = glob.glob('%s/te_containing_fq/*.read_repeat_name.txt' %('/'.join(top_dir)))
    read_repeat = read_repeat_name(read_repeat_files)

    ##cluster reads around insertions
    find_insertion_cluster_bam(align_file, read_repeat, usr_target, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_intact, existingTE_found, teSupportingReads)

    ##output absence
    write_output(top_dir, result, read_repeat, usr_target, exper, TE, required_reads, required_left_reads, required_right_reads, teInsertions, teInsertions_reads, teSupportingReads, existingTE_intact, existingTE_found, teReadClusters, bedtools, lib_size)

 
if __name__ == '__main__':
    main()

