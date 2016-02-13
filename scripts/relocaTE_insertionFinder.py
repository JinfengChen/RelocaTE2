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
    ofile = open(outfile, 'w')
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
                te_name= unit[0]
                m = r_pos.search(unit[4])
                if m:
                    start = m.groups(0)[0]
                    end   = m.groups(0)[1]
                r_count = re.sub(r'\D+', '', unit[7])
                l_count = re.sub(r'\D+', '', unit[8])
                r_supp  = re.sub(r'\D+', '', unit[10])
                l_supp  = re.sub(r'\D+', '', unit[11])
                r_id    = 'repeat_%s_%s_%s' %(chro, start, end)
                print >> ofile, '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s;TSD=%s;Note=%s;Right_junction_reads=%s;Left_junction_reads=%s;Right_support_reads=%s;Left_support_reads=%s;' %(chro, 'RelocaTE2', unit[2], start, end, strand, r_id, te_name, unit[1], ins_type, r_count, l_count, r_supp, l_supp)
    ofile.close()

#[name, seq, start, strand]
def Supporting_count(event, tsd_start, teSupportingReads):
    total = 0
    right = 0
    left  = 0
    read1 = []
    read2 = []
    if teSupportingReads.has_key(event):
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
        #can not delete here, delete after every event finish
        #del teSupportingReads[event]
        return (total, left, right, ','.join(read1), ','.join(read2))
    else:
        return (0,0,0,'','') 

def read_repeat_name(infiles):
    data = defaultdict(list)
    for infile in infiles:
        #print 'read repeat name file: %s' %(infile)
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

def write_output(top_dir, result, read_repeat, usr_target, exper, TE, required_reads, required_left_reads, required_right_reads, teInsertions, teInsertions_reads, teSupportingReads, existingTE_inf, teReadClusters, bedtools, lib_size, teLowQualityReads, teFullReads):
    createdir(result)
    ###remove insertion 
    existingTE_bed     = '%s/existingTE.bed' %('/'.join(top_dir))
    existingTE_bed_chr = '%s/existingTE.%s.bed' %('/'.join(top_dir), usr_target)
    nonref     = '%s/%s.%s.all_nonref_insert.txt' %(result, usr_target, TE)
    nonref_gff = '%s/%s.%s.all_nonref_insert.gff' %(result, usr_target, TE)
    nonsup     = '%s/%s.%s.all_nonref_supporting.txt' %(result, usr_target, TE)
    nonsup_gff = '%s/%s.%s.all_nonref_supporting.gff' %(result, usr_target, TE)
    NONREF     = open ('%s/%s.%s.all_nonref_insert.txt' %(result, usr_target, TE), 'w')
    NONSUP     = open ('%s/%s.%s.all_nonref_supporting.txt' %(result, usr_target, TE), 'w')
    READS      = open ('%s/%s.%s.reads.list' %(result, usr_target, TE), 'w')
    #teInsertions[event][TSD_seq][TSD_start]['count']   += 1   ## total junction reads
    #teInsertions[event][TSD_seq][TSD_start][pos]       += 1   ## right/left junction reads
    #teInsertions[event][TSD_seq][TSD_start][TE_orient] += 1   ## plus/reverse insertions 
    for event in sorted(teInsertions.keys(), key=int):
        print 'event: %s' %(event)
        cluster_collection = []
        start_collection = []
        start_both_junction = 0
        chro = teReadClusters[event]['read_inf']['seq']['chr']
        for start in sorted(teInsertions[event].keys(), key=int):
            #tir1_end = int(start)
            #tir2_end = int(start) - 1
            print 'start: %s' %(start)
            #####supporting reads
            total_supporting, left_supporting, right_supporting = [0, 0, 0]
            left_reads, right_reads = ['', '']
            total_supporting, left_supporting, right_supporting, left_reads, right_reads = Supporting_count(event, start, teSupportingReads)
            repeat_supporting = insertion_family_supporting('%s,%s' %(left_reads, right_reads), read_repeat)
            #del teSupportingReads[event] 
            #####tsd within one start position into one start position
            tsd_count         = {}
            TE_orient, repeat_junction = ['','']
            TE_orient_foward, TE_orient_reverse    = [0,0]
            total_count, left_count, right_count   = [0,0,0]
            total_valid, left_valid, right_valid   = [0,0,0]
            total_low, left_low, right_low = [0,0,0]
            reads             = []
            left_jun_reads    = []
            right_jun_reads = []
            #if two or more TSD at this loci, some might be caused by sequencing error, we add the count
            #and use high frequent one to get te orientation and family
            #if only one TSD found, we acturally use only this one
            for foundTSD in sorted(teInsertions[event][start].keys()):
                total_count += teInsertions[event][start][foundTSD]['count']
                left_count  += teInsertions[event][start][foundTSD]['left']
                right_count += teInsertions[event][start][foundTSD]['right']
                #total_low   += teInsertions[event][start][foundTSD]['count_low']
                #left_low    += teInsertions[event][start][foundTSD]['left_low']
                #right_low   += teInsertions[event][start][foundTSD]['right_low']
                TE_orient_foward  += teInsertions[event][start][foundTSD]['+']
                TE_orient_reverse += teInsertions[event][start][foundTSD]['-']
                tsd_count[foundTSD] = teInsertions[event][start][foundTSD]['count']
                reads.extend(teInsertions_reads[event][start][foundTSD]['read'])
                left_jun_reads.extend(teInsertions_reads[event][start][foundTSD]['left_read'])
                right_jun_reads.extend(teInsertions_reads[event][start][foundTSD]['right_read'])
                print 'TSD count: %s, %s' %(foundTSD, tsd_count[foundTSD])

            tsd_top = OrderedDict(sorted(tsd_count.items(), key=lambda x: x[1])).keys()[-1]
            foundTSD= tsd_top 
            TE_orient         = '+' if int(TE_orient_foward) > int(TE_orient_reverse) else '-'
            repeat_junction   = insertion_family(','.join(reads), read_repeat)

            print 'tsd_top: %s, %s' %(tsd_top, tsd_count[tsd_top])
            print 'reads: %s' %(','.join(reads))
            print 'left: %s' %(','.join(left_jun_reads))
            print 'right: %s' %(','.join(right_jun_reads))
            print 'left_s: %s' %(','.join(left_reads))
            print 'right_s: %s' %(','.join(right_reads))

            remove_start   = 0
            fullreads_flag = 0
            #filter out junctions that fullreads are properly mapped on reference
            left_f, right_f, left_t, right_t = junction_full_reads(left_jun_reads, right_jun_reads, teFullReads)
            print 'fullreads left/right, total reads left/right: %s\t%s\t%s\t%s' %(left_f, right_f, left_t, right_t)
            if float(left_f) >= 0.3 * float(left_t) and float(right_f) >= 0.3 * float(right_t):
                ##junction with reads that properly mapped to reference, should be removed
                #del teInsertions[event][start]
                print 'fullreads remove'
                fullreads_flag = 1

            #filter out these start with junctions only suported by low quality reads
            #real number of high quality reads from sub_cluster not whole cluster
            lowquality_flag = 0
            total_real, total_valid, left_valid, right_valid = lowquality_reads(left_jun_reads, right_jun_reads, teLowQualityReads)
            print '%s\t%s\t%s\t%s' %(total_real, total_valid, left_valid, right_valid)
            #total_valid = total_count - total_low
            #left_valid  = left_count  - left_low
            #right_valid = right_count - right_low
            if total_real == 1 and total_valid == 0:
                ##singleton junction and read is low quality, remove
                print "singleton junction and read is low quality, remove"
                #del teInsertions[event][start]
                lowquality_flag = 1
            elif left_valid == 0 and right_valid == 0:
                ##junction only supported by low quality read, remove
                print "junction only supported by low quality read, remove"
                #del teInsertions[event][start]
                lowquality_flag = 1
            if (fullreads_flag == 0 and lowquality_flag == 0):
            #if (fullreads_flag == 0 and lowquality_flag == 0) or (foundTSD != 'UNK' and lowquality_flag == 0):
                #####information on each start position
                print "store information on each start position"
                start_both_junction = start_both_junction + 1 if (int(left_count) > 0 and int(right_count) > 0) else start_both_junction
                start_collection.append([start, foundTSD, total_count, left_count, right_count, repeat_junction, ','.join(reads), total_supporting, left_supporting, right_supporting, left_reads, right_reads, repeat_supporting, TE_orient])
            else:
                del teInsertions[event][start]
 
        if teSupportingReads.has_key(event):
            del teSupportingReads[event]

        if not len(start_collection) > 0:
            #skip this event if there is no valid junctions
            print "skip this event if there is no valid junctions"
            continue

        if len(teInsertions[event].keys()) > 1:
            #two or more start found in one read cluster
            if int(start_both_junction) == len(teInsertions[event].keys()):
                #we keep all if they all supported by junction at both direction
                print 'we keep all if they all supported by junction at both direction'
                cluster_collection.extend(start_collection)
            else:
                #not all supported at both end of junction
                print 'not all supported at both end of junction'
                if int(start_both_junction) > 0:
                    #we have tsd supported by junction at both direction
                    #we only keep these supported by junction at both direction
                    #even some supported by more than 1 junction reads, we thought due to sequence context
                    #these insertions should have junction at both direction if they are true
                    #now include these single junctions with more than 3 junctions reads support becasue we have more filters before and after
                    print 'we have tsd supported by junction at both direction'
                    for start_it in start_collection:
                        if int(start_it[3]) > 0 and int(start_it[4]) > 0:
                            #both junction
                            cluster_collection.append(start_it)
                        elif int(start_it[3]) >= 3 or int(start_it[4]) >= 3:
                            #single junction
                            cluster_collection.append(start_it)
                        else:
                            continue
                            #for these cluster, we automaticly remove these with only one junction
                            #due to the local coverage, if we found both junction for one insertion
                            #we should find both junction for others too
                            #only remove these low supporting junctions in this case
                            
                else:
                    #do not have tsd supported by junction at both direction
                    #we keep the top one with higher number of junction reads
                    print 'do not have tsd supported by junction at both direction'
                    top_start_it = []
                    top_start_junction = 0
                    for start_it in start_collection:
                        tir1_end = int(start_it[0]) 
                        tir2_end = int(start_it[0]) - 1
                        if existingTE_inf[chro]['start'].has_key(tir1_end) or existingTE_inf[chro]['end'].has_key(tir2_end):
                            continue
                            #remove if overlap with existingTE
                        else:
                            #if int(top_start_junction) < int(start_it[2]):
                            #    top_start_junction = int(start_it[2])
                            #    top_start_it       = start_it
                            #keep all these single junctions that have more than three supporting junctions
                            if int(start_it[2]) >= 3:
                                cluster_collection.append(start_it)
                    #if top_start_junction > 0:
                    #    cluster_collection.append(top_start_it)
                    #else:
                    #    ##all removed due to overlap with existingTE
                    #    continue
                    #    #del teInsertions[event]
        else:
            #only one tsd found
            print "only one tsd found"
            if int(start_collection[0][3]) > 0 and int(start_collection[0][4]) > 0:
                #both junction found, must be insertion
                print "both junction found, must be insertion"
                cluster_collection.extend(start_collection)
            else:
                #only one junction, check if overlap with existingTE
                print "only one junction, check if overlap with existingTE"
                tir1_end = int(start_collection[0][0])
                tir2_end = int(start_collection[0][0]) - 1
                if existingTE_inf[chro]['start'].has_key(tir1_end) or existingTE_inf[chro]['end'].has_key(tir2_end):
                    #overlap with existingTE, delete record from teInsertions
                    print "overlap with existingTE, delete record from teInsertions"
                    continue
                    #del teInsertions[event]
                else:
                    #not overlap with existingTE, store for output
                    print "not overlap with existingTE, store for output"
                    cluster_collection.extend(start_collection)
         
        #[start, foundTSD, total_count, left_count, right_count, repeat_junction, ','.join(reads),
        #total_supporting, left_supporting, right_supporting, left_reads, right_reads, repeat_supporting, TE_orient]
        for insertion in cluster_collection:
            #####insertion te name
            r_supporting = insertion[12]
            r_junction   = insertion[5]
            repeat_family     = 'NA'
            if r_junction == r_supporting and r_junction != 'NA':
                repeat_family = r_junction
            elif r_junction != r_supporting and r_junction != 'NA' and r_supporting != 'NA':
                repeat_family = '%s/%s' %(r_junction, r_supporting)
            elif r_junction != 'NA':
                repeat_family = r_junction
            elif r_supporting != 'NA':
                repeat_family = r_supporting
            
            #####write to file
            t_count = insertion[2]
            l_count = insertion[3]
            r_count = insertion[4]
            i_start = insertion[0]
            i_tsd   = insertion[1]
            i_reads = insertion[6]
            l_reads = insertion[10]
            r_reads = insertion[11]
            t_supporting = insertion[7]
            l_supporting = insertion[8]
            r_supporting = insertion[9]
            t_orient     = insertion[13]
            
            if int(l_count) >= int(required_left_reads) or int(r_count) >= int(required_right_reads):
                #at lease need one end have junction reads
                print "at lease need one end have junction reads"
                coor       = int(i_start) + (len(i_tsd) - 1)
                coor_start = coor - (len(i_tsd) - 1)
                print >> READS, '%s\t%s:%s..%s\tJunction_reads\t%s' %(TE, usr_target, coor_start, coor, i_reads)
                print >> READS, '%s\t%s:%s..%s\tLeft_supporting_reads\t%s' %(TE, usr_target, coor_start, coor, l_reads)
                print >> READS, '%s\t%s:%s..%s\tRight_supporting_reads\t%s' %(TE, usr_target, coor_start, coor, r_reads)
                if int(l_count) > 0 and int(r_count) >0:
                    #both ends have junction reads
                    print "both ends have junction reads"
                    print >> NONREF, '%s\t%s\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, i_tsd, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                else:
                    #only have end with junction
                    if int(r_supporting) >= 1 and int(l_supporting) >= 1:
                        #only one end with junction but both end with supporting reads
                        print "only one end with junction but both end with supporting reads"
                        print >> NONREF, '%s\tsupporting_junction\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                    else:
                        #only one end with junction and only one end with supporting reads
                        print "only one end with junction and only one end with supporting reads"
                        if (int(r_supporting) >= 1 and int(l_count) >= 1) or (int(l_supporting) >= 1 and int(r_count) >= 1):
                            print "supporting junction"
                            print >> NONREF, '%s\tsupporting_junction\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                        elif int(t_supporting) + int(t_count) == 1:
                            #singleton
                            print "singleton"
                            print >> NONREF, '%s\tsingleton\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                        else:
                            #insufficient
                            print "insufficient"
                            print >> NONREF, '%s\tinsufficient_data\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
            else:
                #no junction reads
                print "no junction reads"
                if int(r_supporting) >= 1 and int(l_supporting) >= 1:
                    #no junction reads, but both end with supporting reads
                    print >> NONREF, '%s\tsupporting_reads\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                elif int(t_supporting) == 1:
                    print >> NONREF, '%s\tsingleton\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                else:
                    print >> NONREF, '%s\tinsufficient_data\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(repeat_family, exper, usr_target, coor_start, coor, t_orient, t_count, r_count, l_count, t_supporting, r_supporting, l_supporting)
                #pass
    ###only supporting reads, no junction record
    for event in sorted(teSupportingReads.keys()):
        if teInsertions.has_key(event):
            #skip any event that have junctions, already output in previous step
            #print >> NONSUP, 'Found' %(str(event))
            continue
        #print >> NONSUP, '>%s' %(str(event))
        sub_events = defaultdict(lambda : list())
        sub_strand = defaultdict(lambda : list())
        sub_side   = defaultdict(lambda : int()) ## record if have read from both strand
        sub_name   = 'repeat_name'
        for read in teSupportingReads[event]:
            name   = read[0]
            seq    = read[1]
            start  = read[2]
            strand = read[3]
            sub_side[strand] = 1
            #print >> NONSUP, '%s\t%s\t%s\t%s' %(name, seq, start, strand)
            #repeat_name, repeat_side = read_to_repeat(name, read_repeat)
            #print >> NONSUP, repeat_name, repeat_side
            #repeat_name = read_repeat[name][1] if read_repeat.has_key(name) else 'NA'
            #repeat_side = read_repeat[name][2] if read_repeat.has_key(name) else 'NA'
            sub_events[strand].append(read)
        if len(sub_events.keys()) == 2:
            ##supporting from both end
            #print >> NONSUP, 'IN'
            
            ins_start = get_boundary(sub_events['+'], 'left')
            l_support = len(sub_events['+'])
            ins_end   = get_boundary(sub_events['-'], 'right')
            r_support = len(sub_events['-'])
            t_support = l_support + r_support
            if ins_start > ins_end:
                continue
            print >> NONSUP, '%s\tsupporting_reads\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(sub_name, exper, usr_target, ins_start, ins_end, '+', '0','0','0', t_support, r_support, l_support)
        elif sub_events.keys()[0] == '+':
            
            ins_start = get_boundary(sub_events['+'], 'left')
            ins_end   = int(ins_start + float(lib_size)*(1 + 0.2)) # insertion size of library * (1 + sd of library)
            l_support = len(sub_events['+'])
            r_support = 0
            t_support = l_support + r_support
            print >> NONSUP, '%s\tsupporting_reads\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(sub_name, exper, usr_target, ins_start, ins_end, '+', '0','0','0', t_support, r_support, l_support)

        elif sub_events.keys()[0] == '-':
       
            ins_end   = get_boundary(sub_events['-'], 'right')
            ins_start = int(ins_end - float(lib_size)*(1 + 0.2)) # insertion size of library * (1 + sd of library)
            l_support = 0
            r_support = len(sub_events['-'])
            t_support = l_support + r_support
            print >> NONSUP, '%s\tsupporting_reads\t%s\t%s\t%s..%s\t%s\tT:%s\tR:%s\tL:%s\tST:%s\tSR:%s\tSL:%s' %(sub_name, exper, usr_target, ins_start, ins_end, '+', '0','0','0', t_support, r_support, l_support)
        else:
            #print >> NONSUP, 'More than one repeatname or 0'
            pass
    NONREF.close()
    NONSUP.close()
    READS.close()
    txt2gff(nonref, nonref_gff, 'Non-reference, not found in reference')
    txt2gff(nonsup, nonsup_gff, 'Non-reference, not found in reference')
    if os.path.isfile(nonsup_gff) and os.path.getsize(nonsup_gff) > 0:
        os.system('grep -P \"%s\\t\" %s > %s' %(usr_target, existingTE_bed, existingTE_bed_chr))
        os.system('%s intersect -v -a %s -b %s >> %s' %(bedtools, nonsup_gff, existingTE_bed_chr, nonref_gff))
        os.system('rm %s' %(existingTE_bed_chr))

def lowquality_reads(left_jun_reads, right_jun_reads, teLowQualityReads):
    real_c, real_t, real_left_t, real_right_t = [0,0,0,0]
    for rd1 in left_jun_reads:
        print 'left reads: %s' %(rd1)
        real_c += 1
        if not teLowQualityReads.has_key(rd1):
            real_left_t   += 1
            real_t        += 1
            print 'left: %s\t%s\t%s' %(real_c, real_left_t, real_t)
    for rd2 in right_jun_reads:
        print 'right reads: %s' %(rd2)
        real_c += 1
        if not teLowQualityReads.has_key(rd2):
            real_right_t  += 1
            real_t        += 1
            print 'right: %s\t%s\t%s' %(real_c, real_right_t, real_t)
    return [real_c, real_t, real_left_t, real_right_t]

def junction_full_reads(left_jun_reads, right_jun_reads, teFullReads):
    left_f, right_f, left_t, right_t = [0,0,0,0]
    for rd1 in left_jun_reads:
        left_t += 1
        if teFullReads.has_key(rd1):
            print 'fullreads: %s' %(rd1)
            left_f   += 1
    for rd2 in right_jun_reads:
        right_t += 1
        if teFullReads.has_key(rd2):
            print 'fullreads: %s' %(rd2)
            right_f  += 1
    return [left_f, right_f, left_t, right_t]


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

def TSD_from_read_depth(r, read_repeat, teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teLowQualityReads):
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
        print 'tsdfiner%s' %(cluster)
        for name in teReadClusters[cluster]['read_inf'].keys():
            if name == 'seq':
                #skip empty line
                continue
            seq    = teReadClusters[cluster]['read_inf'][name]['seq']
            start  = teReadClusters[cluster]['read_inf'][name]['start']
            end    = teReadClusters[cluster]['read_inf'][name]['end']
            strand = teReadClusters[cluster]['read_inf'][name]['strand']
            print name, start, end, seq, strand
            #end    = int(start) + len(seq)
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
        print len(left_reads.keys()), len(right_reads.keys())
        if (len(left_reads.keys()) > 1 and len(right_reads.keys()) >= 1) or (len(left_reads.keys()) >= 1 and len(right_reads.keys()) > 1):
            ##more than two junction
            print 'more than two junction'
            count_tsd = 0
            pairs_tsd = defaultdict(lambda : int())
            ##find pairs for one insertions
            print 'find pairs for one insertions'
            for start1 in sorted(left_reads.keys(), key=int):
                min_dist = -100
                min_pair = ''
                for start2 in sorted(right_reads.keys(), key=int):
                    if min_dist < 0:
                        min_dist = abs(int(start2) - int(start1))
                        min_pair = start2
                    elif min_dist > abs(int(start2) - int(start1)):
                        #20160212 modifed to deal with both close insertions and false mapping reads
                        if min_dist <= 100 and int(start2) < int(start1):
                            #already have one pair of TSD
                            if len(right_reads[start2]) >= 2 and len(right_reads[min_pair]) >= 2:
                                #both right_reads have more than 2 reads support, likely true. we pick the one with short distance
                                min_dist = abs(int(start2) - int(start1))
                                min_pair = start2
                            else:
                                #not all have more than 2 reads support, likely one is false. we pick the one with higher depth
                                if len(right_reads[start2]) > len(right_reads[min_pair]):
                                    min_dist = abs(int(start2) - int(start1))
                                    min_pair = start2
                        else:
                            min_dist = abs(int(start2) - int(start1))
                            min_pair = start2
                        #20160212
                    #elif min_dist > abs(int(start2) - int(start1)):
                    #    if min_dist <= 100 and int(start2) < int(start1):
                    #        ##likely start1 and start2 overlap (find TSD), we choose one with more sequence depth
                    #        if len(right_reads[start2]) > len(right_reads[min_pair]):
                    #            min_dist = abs(int(start2) - int(start1))
                    #            min_pair = start2
                    #    else:
                    #        ##may not TSD
                    #        min_dist = abs(int(start2) - int(start1))
                    #        min_pair = start2
                if min_dist <= 100:
                    ##find pairs
                    print 'find pairs'
                    count_tsd += 1
                    pairs_tsd[start1]   = 1
                    pairs_tsd[min_pair] = 1
                    teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 2
                    for read in left_reads[start1]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                    for read in right_reads[min_pair]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                else:
                    ##do not find pairs
                    print 'do not find pairs'
                    count_tsd += 1
                    pairs_tsd[start1]   = 1
                    teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                    for read in left_reads[start1]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
            ##set unpaired
            print 'set unpaired'
            for start2 in right_reads.keys():
                if not pairs_tsd.has_key(start2):
                    #not paired junction
                    count_tsd += 1
                    teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                    for read in right_reads[start2]:
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                        teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                        calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
        elif len(left_reads.keys()) > 1:
            #more than two left junctions
            print 'more than two left junctions'
            count_tsd = 0
            for start1 in left_reads.keys():
                count_tsd += 1
                teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
        elif len(right_reads.keys()) > 1:
            #more than two right junction
            print 'more than two right junction'
            count_tsd = 0
            for start2 in right_reads.keys():
                count_tsd += 1
                teReadClusters_sub_type['%s-%s' %(cluster, count_tsd)] = 1
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                    teReadClusters_sub['%s-%s' %(cluster, count_tsd)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-%s' %(cluster, count_tsd), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
        elif len(left_reads.keys()) == 1 and len(right_reads.keys()) == 1:
            ##one right and one left junction
            print 'one right and one left junction'
            print left_reads.keys()[0], right_reads.keys()[0]
            if abs(int(left_reads.keys()[0]) - int(right_reads.keys()[0])) > 100:
                ##two junctions are far from each other, might be one end from two insertion
                print 'two junctions are far from each other, might be one end from two insertion'
                teReadClusters_sub_type['%s-1' %(cluster)] = 1
                teReadClusters_sub_type['%s-2' %(cluster)] = 1
                start1 = left_reads.keys()[0]
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                start2 = right_reads.keys()[0]
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-2' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-2' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-2' %(cluster)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                    teReadClusters_sub['%s-2' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-2' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
            else:
                ##one junction
                print 'two junctions are for one insertion'
                teReadClusters_sub_type['%s-1' %(cluster)] = 2
                start1 = left_reads.keys()[0]
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
                start2 = right_reads.keys()[0]
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth) 
        else:
            ##one junction with one end support
            print 'one junction with one end support'
            teReadClusters_sub_type['%s-1' %(cluster)] = 1
            if len(left_reads.keys()) > 0:
                #teReadClusters_sub_type['%s-1' %(cluster)] = 1
                start1 = left_reads.keys()[0] 
                for read in left_reads[start1]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)
            elif len(right_reads.keys()) > 0:
                start2 = right_reads.keys()[0]
                for read in right_reads[start2]:
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['seq']    = teReadClusters[cluster]['read_inf'][read]['seq']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['start']  = teReadClusters[cluster]['read_inf'][read]['start']
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['end']  = teReadClusters[cluster]['read_inf'][read]['end']  
                    teReadClusters_sub['%s-1' %(cluster)]['read_inf'][read]['strand'] = teReadClusters[cluster]['read_inf'][read]['strand']
                    calculate_cluster_depth('%s-1' %(cluster), teReadClusters[cluster]['read_inf'][read]['seq'], teReadClusters[cluster]['read_inf'][read]['start'], teReadClusters[cluster]['read_inf'][read]['end'], read, teReadClusters[cluster]['read_inf'][read]['strand'], teReadClusters_sub, teReadClusters_sub_count, teReadClusters_sub_depth)

        #print 'TSD finder: %s' %(cluster)
        ###teReadCluster_sub_depth add above
        ###deal with teReadClusters_sub, still store at cluster in teInsertions, write a TSD_check for this only.
        for sub_cluster in teReadClusters_sub_depth.keys():
            print 'sub_cluster: %s' %(sub_cluster)
            
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
                    end    = teReadClusters_sub[sub_cluster]['read_inf'][name]['end']
                    strand = teReadClusters_sub[sub_cluster]['read_inf'][name]['strand']
                    TSD_check_single(cluster, seq, chro, start, end, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teLowQualityReads)
            else:
                #estimate TSD_len using read depth
                TSD_len = 0
                if tsd_finder(sub_cluster, 1, teReadClusters_sub_count, teReadClusters_sub_depth):
                    TSD_len = tsd_finder(sub_cluster, 1, teReadClusters_sub_count, teReadClusters_sub_depth)
                elif tsd_finder(sub_cluster, 0.8, teReadClusters_sub_count, teReadClusters_sub_depth):
                    TSD_len = tsd_finder(sub_cluster, 0.8, teReadClusters_sub_count, teReadClusters_sub_depth)
                elif tsd_finder(sub_cluster, 0.6, teReadClusters_sub_count, teReadClusters_sub_depth):
                    TSD_len = tsd_finder(sub_cluster, 0.6, teReadClusters_sub_count, teReadClusters_sub_depth)

                #estimate TSD_len using position of left and right read cluster
                TSD_len_1 = 0
                TSD_len_1 = TSD_len_calculate(teReadClusters_sub, sub_cluster, r)
                

                #merge, use TSD_len_1 if not consistent
                print 'TSD_merge: TSD_len, TSD_len_1 %s\t%s' %(TSD_len, TSD_len_1)
                if not int(TSD_len) == int(TSD_len_1):
                    TSD_len = int(TSD_len_1)               

                if TSD_len > 0:
                    print 'TSD found: %s' %(TSD_len)
                    TSD = '.'*TSD_len
                    for name1 in teReadClusters_sub[sub_cluster]['read_inf'].keys():
                        real_name = r.search(name1).groups(0)[0] if r.search(name1) else ''
                        seq    = teReadClusters_sub[sub_cluster]['read_inf'][name1]['seq']
                        start  = teReadClusters_sub[sub_cluster]['read_inf'][name1]['start']
                        end  = teReadClusters_sub[sub_cluster]['read_inf'][name1]['end']
                        strand = teReadClusters_sub[sub_cluster]['read_inf'][name1]['strand']
                        #print name1, seq, start, strand, chro
                        TSD_check_cluster(cluster, seq, chro, start, end, real_name, read_repeat, name1, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teLowQualityReads)
                        #print 'after tsd_check'
                else:          
                    #what if we can not find TSD? still could be insertions
                    TSD = 'UKN'
                    print 'no TSD found: %s' %(TSD)
                    for name in teReadClusters_sub[sub_cluster]['read_inf'].keys():
                        real_name = r.search(name).groups(0)[0] if r.search(name) else ''
                        seq    = teReadClusters_sub[sub_cluster]['read_inf'][name]['seq']
                        start  = teReadClusters_sub[sub_cluster]['read_inf'][name]['start']
                        end    = teReadClusters_sub[sub_cluster]['read_inf'][name]['end']
                        strand = teReadClusters_sub[sub_cluster]['read_inf'][name]['strand']
                        TSD_check_cluster(cluster, seq, chro, start, end, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teLowQualityReads)
            #print 'End of cycle'

def tsd_finder(sub_cluster, tsd_depth, teReadClusters_sub_count, teReadClusters_sub_depth):
    TSD_len  = 0
    read_total = teReadClusters_sub_count[sub_cluster]['read_count']
    for chrs_pos in sorted(teReadClusters_sub_depth[sub_cluster]['read_inf']['depth'].keys(), key=int):
        depth = teReadClusters_sub_depth[sub_cluster]['read_inf']['depth'][chrs_pos]
        if float(depth) >= float(tsd_depth)*float(read_total):
            TSD_len += 1
    return TSD_len

def TSD_len_calculate(teReadClusters_sub, sub_cluster, r):
    ##seq is entire trimmd read, not just the TSD portion of the read
    ##start is the first postition of the entire read match to ref
    TSD_len_1 = 0
    r5 = re.compile(r'start:[53]$')
    r3 = re.compile(r'end:[53]$')
    #left and right end of TSD on chromosome
    TSD_left  = defaultdict(lambda : int())
    TSD_right = defaultdict(lambda : int())
    #print 'TSD_check\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    ##start means that the TE was removed from the start of the read
    ##5 means the trimmed end mapps to the 5prime end of the TE
    ##3 means the trimmed end mapps to the 3prime end of the TE
    print 'TSD_len_calculate: %s' %(sub_cluster)
    for name in teReadClusters_sub[sub_cluster]['read_inf'].keys():
        real_name = r.search(name).groups(0)[0] if r.search(name) else ''
        seq    = teReadClusters_sub[sub_cluster]['read_inf'][name]['seq']
        start  = teReadClusters_sub[sub_cluster]['read_inf'][name]['start']
        end    = teReadClusters_sub[sub_cluster]['read_inf'][name]['end']
        strand = teReadClusters_sub[sub_cluster]['read_inf'][name]['strand']
        pos       = ''
        TE_orient = 0
        if strand == '+':
            if r5.search(name):
                pos       = 'right'
                TE_orient = '-' if name[-1] == '5' else '+'
                TSD_left[int(start)] += 1
            elif r3.search(name):
                pos       = 'left'
                TE_orient = '+' if name[-1] == '5' else '-'
                #TSD_right[int(start) + len(seq) - 1] += 1
                TSD_right[int(end) - 1]  += 1
        elif strand == '-':
            if r5.search(name):
                pos       = 'left'
                TE_orient = '+' if name[-1] == '5' else '-'
                #TSD_right[int(start) + len(seq) - 1] += 1
                TSD_right[int(end) - 1]  += 1
            elif r3.search(name):
                pos       = 'right'
                TE_orient = '-' if name[-1] == '5' else '+'
                TSD_left[int(start)] += 1
        if pos == 'right':
            print '%s\t%s\t%s\t%s\t%s' %(name, start, pos, len(TSD_left.keys()), len(TSD_right.keys()))
        elif pos == 'left':
            print '%s\t%s\t%s\t%s\t%s' %(name, int(start) + len(seq) - 1, pos, len(TSD_left.keys()), len(TSD_right.keys()))
            print '%s\t%s\t%s\t%s\t%s' %(name, int(end), pos, len(TSD_left.keys()), len(TSD_right.keys()))
    TSD_left_sort  = OrderedDict(sorted(TSD_left.items(), key=lambda x: x[1])) 
    TSD_right_sort = OrderedDict(sorted(TSD_right.items(), key=lambda x: x[1]))
    print 'sorted max value: %s\t%s' %(int(TSD_left_sort.values()[-1]), int(TSD_right_sort.values()[-1]))
    print 'sorted max key: %s\t%s' %(int(TSD_left_sort.keys()[-1]), int(TSD_right_sort.keys()[-1]))
    #more than two read support the TSD boundary
    if int(TSD_left_sort.values()[-1]) >= 1 and int(TSD_right_sort.values()[-1]) >= 1:
        TSD_len_1 = int(TSD_right_sort.keys()[-1]) - int(TSD_left_sort.keys()[-1]) + 1
    else:
        TSD_len_1 = 0
    print 'Estimated TSD_len_1: %s' %(TSD_len_1)
    return TSD_len_1

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
                calculate_cluster_depth(count, seq, start, end, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
        elif not r.search(name) and not record.is_paired:
            #reads not matched to repeat and not mates of junctions
            #reads are mates of reads matched to middle of repeat
            #supporting reads
            teSupportingReads[count].append([name, seq, start, strand])
        #elif not r.search(name) and record.is_paired:
            #reads not matched to repeat but mates of junctions
        #    teSupportingReads[count].append([name, seq, start, strand])
    else:
        #if start and end do not fall within last start and end
        #we now have a different insertion event
        count += 1
        if r.search(name):
            real_name = r.search(name).groups(0)[0]
            if not r_tsd.search(TSD):
                TSD_check(count, seq, chro, start, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
            else:
                calculate_cluster_depth(count, seq, start, end, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
        elif not r.search(name) and not record.is_paired:
            #reads not matched to repeat and not mates of junctions
            #reads are mates of reads matched to middle of repeat
            #supporting reads
            teSupportingReads[count].append([name, seq, start, strand])
        #elif not r.search(name) and record.is_paired:
            #reads not matched to repeat but mates of junctions
        #    teSupportingReads[count].append([name, seq, start, strand])
        #initial insertion site boundary
        bin_ins = [int(start), int(end)]
        #print '%s\t%s' %(count, bin_ins)
    return (bin_ins, count)

def existingTE(infile, existingTE_inf, existingTE_found):
    r = re.compile(r'(\S+)\t(\S+):(\d+)\.\.(\d+)')
    if infile != 'NONE':
        blat = 0
        with open(infile, 'r') as filehd:
            for line in filehd:
                line = line.rstrip() 
                if not blat:
                    if r.search(line):
                        te, chrs, start, end = r.search(line).groups(0)
                        start, end = sorted([start, end], key=int)
                        existingTE_inf[te]['start'][int(start)] = line
                        existingTE_inf[te]['end'][int(end)]     = line
                        existingTE_found[line]['start']= 0
                        existingTE_found[line]['end']  = 0


def existingTE_RM_ALL(top_dir, infile, existingTE_inf, chro):
    #ofile_RM = open('%s/existingTE.bed' %(top_dir), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\s+',line)
                if not unit[0] == '':
                    unit.insert(0, '')
                #print line
                #print unit[6], unit[7], unit[9], unit[12], unit[13], unit[14]
                if not unit[5] == chro:
                    continue
                if unit[9] == '+':
                    for i in range(int(unit[6])-2, int(unit[6])+3):
                        existingTE_inf[unit[5]]['start'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7], '1', '+')
                        #print unit[10], 'start', unit[6]
                    for i in range(int(unit[7])-2, int(unit[7])+3):
                        existingTE_inf[unit[5]]['end'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '+')
                        #print unit[10], 'end', unit[7]

                    ##if this repeat is a intact element
                    #intact = 0
                    #if int(unit[12]) == 1 and len(unit[14]) == 3:
                    #    unit[14] =re.sub(r'\(|\)', '', unit[14])
                    #    if int(unit[14]) == 0:
                    #        intact = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])), str(int(unit[7])), unit[11],unit[6],unit[7], intact, '+')
                elif unit[9] == 'C':
                    for i in range(int(unit[6])-2, int(unit[6])+3):
                        existingTE_inf[unit[5]]['start'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])-2), str(int(unit[6])+2), unit[11],unit[6],unit[7],'1', '-')
                        #print unit[10], 'start', unit[6]
                    for i in range(int(unit[7])-2, int(unit[7])+3):
                        existingTE_inf[unit[5]]['end'][int(i)] = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[7])-2), str(int(unit[7])+2), unit[11],unit[6],unit[7], '1', '-')
                        #print unit[10], 'end', unit[7]
                    #intact = 0
                    #if int(unit[14]) == 1 and len(unit[12]) == 3:
                    #    unit[12] =re.sub(r'\(|\)', '', unit[12])
                    #    if int(unit[12]) == 0:
                    #        intact = 1
                    #print >> ofile_RM, '%s\t%s\t%s\t%s:%s-%s\t%s\t%s' %(unit[5], str(int(unit[6])), str(int(unit[7])), unit[11],unit[6],unit[7], intact, '-')
    #ofile_RM.close()

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

def calculate_cluster_depth(event, seq, start, end, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth):
    teReadClusters_count[event]['read_count'] += 1
    teReadClusters[event]['read_inf'][name]['seq']   = seq
    teReadClusters[event]['read_inf'][name]['start'] = start
    teReadClusters[event]['read_inf'][name]['end'] = end
    teReadClusters[event]['read_inf'][name]['strand']= strand
    for i in range(int(start), int(end)):
        teReadClusters_depth[event]['read_inf']['depth'][i] += 1

def TSD_check_single(event, seq, chro, start, end, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teLowQualityReads):
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
    if verbose > 2: print 'TSD_check_single\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
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
            #TSD_start = int(start) + len(seq)
            TSD_start = int(end)
    elif strand == '-':
        if r5.search(name):
            result    = 1
            TSD_seq   = 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            #TSD_start = int(start) + len(seq)
            TSD_start = int(end)
        elif r3.search(name):
            result    = 1
            TSD_seq   = 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
    if verbose > 2: print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    if result and TE_orient:
        tir1_end, tir2_end = [0, 0]
        #if 0:
        #    continue
        if pos == 'left':
            tir1_end = int(TSD_start)
            if verbose > 2: print 'tir1: %s' %(tir1_end)
        elif pos == 'right':
            tir2_end = int(TSD_start) - 1
            if verbose > 2: print 'tir2: %s' %(tir2_end)
        if tir1_end > 0 and existingTE_inf[chro]['start'].has_key(tir1_end):
            te_id = existingTE_inf[chro]['start'][tir1_end]
        #    existingTE_found[te_id]['start'] += 1
            if verbose > 2: print 'tir1'
        elif tir2_end > 0 and existingTE_inf[chro]['end'].has_key(tir2_end):
            te_id = existingTE_inf[chro]['end'][tir2_end]
        #    existingTE_found[te_id]['end'] += 1
            if verbose > 2: print 'tir2'
        else:
            if verbose > 2: print 'not match'
            ##non reference insertions
            if teLowQualityReads.has_key(name):
                teInsertions[event][TSD_start][TSD_seq]['count_low']     += 1
                teInsertions[event][TSD_start][TSD_seq]['%s_low' %(pos)] += 1            

            teInsertions[event][TSD_start][TSD_seq]['count']   += 1   ## total junction reads
            teInsertions[event][TSD_start][TSD_seq][pos]       += 1   ## right/left junction reads
            teInsertions[event][TSD_start][TSD_seq][TE_orient] += 1   ## plus/reverse insertions
            #read_name = re.sub(r':start|:end', '', name)
     
            if pos == 'left':
                teInsertions_reads[event][TSD_start][TSD_seq]['left_read'].append(name)
            elif pos == 'right':
                teInsertions_reads[event][TSD_start][TSD_seq]['right_read'].append(name)

            teInsertions_reads[event][TSD_start][TSD_seq]['read'].append(name)
            #print '1: %s\t 2: %s' %(read_name, teInsertions_reads[event][TSD_seq][TSD_start]['read'])
            if verbose > 2: print 'C: %s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient)

def TSD_check(event, seq, chro, start, end, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found):
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
    if verbose > 2: print 'TSD_check\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
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
            #TSD_start = int(start) + (len(seq)-len(TSD))
            TSD_start = int(end)  - len(TSD)
    elif strand == '-':
        if r5.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
            result    = 1
            TSD_seq   = r3_tsd.search(seq).groups(0)[0] if r3_tsd.search(seq) else 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            #TSD_start = int(start) + (len(seq)-len(TSD))
            TSD_start = int(end)  - len(TSD)
        elif r3.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
            result    = 1
            TSD_seq   = r5_tsd.search(seq).groups(0)[0] if r5_tsd.search(seq) else 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
    if verbose > 2: print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    if result and TE_orient:
        tir1_end, tir2_end = [0, 0]
        if 0:
            continue
        if pos == 'left':
            #tir1_end = int(start) + len(seq)
            tir1_end = int(end)
            if verbose > 2: print 'tir1: %s' %(tir1_end)
        elif pos == 'right':
            tir2_end = int(start) - 1
            if verbose > 2: print 'tir2: %s' %(tir2_end)
        if tir1_end > 0 and existingTE_inf[chro]['start'].has_key(tir1_end):
            te_id = existingTE_inf[chro]['start'][tir1_end]
        #    #existingTE_found[te_id]['start'] += 1
            if verbose > 2: print 'tir1'
        elif tir2_end > 0 and existingTE_inf[chro]['end'].has_key(tir2_end):
            te_id = existingTE_inf[chro]['end'][tir2_end]
        #    #existingTE_found[te_id]['end'] += 1
            if verbose > 2: print 'tir2'
        else:
            #print 'not match'
            ##non reference insertions
            teInsertions[event][TSD_start][TSD_seq]['count']   += 1   ## total junction reads
            teInsertions[event][TSD_start][TSD_seq][pos]       += 1   ## right/left junction reads
            teInsertions[event][TSD_start][TSD_seq][TE_orient] += 1   ## plus/reverse insertions
            #read_name = re.sub(r':start|:end', '', name)
            teInsertions_reads[event][TSD_start][TSD_seq]['read'].append(name)
            #print '1: %s\t 2: %s' %(read_name, teInsertions_reads[event][TSD_seq][TSD_start]['read'])
            if verbose > 2: print 'C: %s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient)

def convert_tag(tag):
    tags = {}
    for t in tag:
        tags[t[0]] = t[1]
    return tags

def TSD_check_cluster(event, seq, chro, start, end, real_name, read_repeat, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teLowQualityReads):
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
    if verbose > 2: print 'TSD_check\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
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
            #TSD_start = int(start) + (len(seq)-len(TSD))
            TSD_start = int(end)  - len(TSD)
    elif strand == '-':
        if r5.search(name) and (r5_tsd.search(rev_com) or r3_tsd.search(seq)):
            result    = 1
            TSD_seq   = r3_tsd.search(seq).groups(0)[0] if r3_tsd.search(seq) else 'UNK'
            pos       = 'left'
            TE_orient = '+' if name[-1] == '5' else '-'
            #TSD_start = int(start) + (len(seq)-len(TSD))
            TSD_start = int(end)  - len(TSD)
        elif r3.search(name) and (r5_tsd.search(seq) or r3_tsd.search(rev_com)):
            result    = 1
            TSD_seq   = r5_tsd.search(seq).groups(0)[0] if r5_tsd.search(seq) else 'UNK'
            pos       = 'right'
            TE_orient = '-' if name[-1] == '5' else '+'
            TSD_start = int(start)
    if verbose > 2: print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient, pos, repeat)
    if result and TE_orient:
        tir1_end, tir2_end = [0, 0]
        if 0: # do not filter exiting TE when insertion have both junction support
            continue
        #if pos == 'left':
        #    tir1_end = int(start) + len(seq)
        #    print 'tir1: %s' %(tir1_end)
        #elif pos == 'right':
        #    tir2_end = int(start) - 1
        #    print 'tir2: %s' %(tir2_end)
        #if tir1_end > 0 and existingTE_inf[chro]['start'].has_key(tir1_end):
        #    te_id = existingTE_inf[chro]['start'][tir1_end]
        #    #existingTE_found[te_id]['start'] += 1
        #    print 'tir1'
        #elif tir2_end > 0 and existingTE_inf[chro]['end'].has_key(tir2_end):
        #    te_id = existingTE_inf[chro]['end'][tir2_end]
        #    #existingTE_found[te_id]['end'] += 1
        #    print 'tir2'
        else:
            #print 'not match'
            ##non reference insertions
            if teLowQualityReads.has_key(name):
                teInsertions[event][TSD_start][TSD_seq]['count_low']     += 1
                teInsertions[event][TSD_start][TSD_seq]['%s_low' %(pos)] += 1            
                
            teInsertions[event][TSD_start][TSD_seq]['count']   += 1   ## total junction reads
            teInsertions[event][TSD_start][TSD_seq][pos]       += 1   ## right/left junction reads
            teInsertions[event][TSD_start][TSD_seq][TE_orient] += 1   ## plus/reverse insertions
            #read_name = re.sub(r':start|:end', '', name)
            if pos == 'left':
                teInsertions_reads[event][TSD_start][TSD_seq]['left_read'].append(name)
            elif pos == 'right':
                teInsertions_reads[event][TSD_start][TSD_seq]['right_read'].append(name)
            teInsertions_reads[event][TSD_start][TSD_seq]['read'].append(name)
            #print '1: %s\t 2: %s' %(read_name, teInsertions_reads[event][TSD_seq][TSD_start]['read'])
            if verbose > 2: print 'C: %s\t%s\t%s\t%s\t%s' %(event, name, TSD_seq, TSD_start, TE_orient)


def find_insertion_cluster_bam(align_file, read_repeat, target, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found, teSupportingReads, teLowQualityReads, teJunctionReads, teFullReads, teReadsInfo, mm_allow):
    r = re.compile(r'(.*):(start|end):(5|3)')
    r_tsd = re.compile(r'UNK|UKN|unknown', re.IGNORECASE)
    r_cg1 = re.compile(r'[S]')
    r_cg2 = re.compile(r'[ID]')
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
            #end    = int(start) + int(length) - 1 #should not allowed for indel or softclip
            end    = int(record.reference_end) + 1
            #print 'end compare %s: %s, %s' %(name, end)
            if verbose > 4: print 'find_insertion_bam: %s\t%s' %(name, start)
            #filter false junctions
            if r.search(name):
                #only for junction reads
                if verbose > 4: print 'Junction: %s\t%s\t%s\t%s' %(name, chro, start, end)
                extend = 10
                read_order = 0
                if record.is_read1:
                    read_order = 1
                elif record.is_read2:
                    read_order = 2
                jun_read_name = r.search(name).groups(0)[0]
                if verbose > 4: print 'Jun_read_name: %s' %(jun_read_name)
                if jun_read_name[-2:] == '/1':
                    if verbose > 4: print 'Jun_read_name: /1'
                    if not teJunctionReads.has_key(jun_read_name[:-2]) or not teJunctionReads[jun_read_name[:-2]].has_key(1):
                        if verbose > 4: print 'Jun_read_name is not a key of teJunctionReads'
                        continue
                    #read has label for read1, we use label for both paired or unpaired reads
                    #if teJunctionReads[jun_read_name[:-2]][1] == 1:
                    ##we assume that the fullread and junction reads mapped to the same position on chromosome, if not we will do nothing
                    if teJunctionReads[jun_read_name[:-2]][1][3] == 1:
                        #fullread mapped properly to genome
                        teFullReads[name] = 1
                    elif teJunctionReads[jun_read_name[:-2]][1][0] == chro and int(teJunctionReads[jun_read_name[:-2]][1][1]) == start:
                        #fullread overlap with junction read
                        if int(teJunctionReads[jun_read_name[:-2]][1][2]) >= end + extend:
                            #fullread extend more than 5 bp into insertion, suggesting no insertion
                            teFullReads[name] = 1
                    elif teJunctionReads[jun_read_name[:-2]][1][0] == chro and int(teJunctionReads[jun_read_name[:-2]][1][2]) == end:
                        if int(teJunctionReads[jun_read_name[:-2]][1][1]) <= start - extend:
                            teFullReads[name] = 1
                    #elif not teJunctionReads[jun_read_name[:-2]][1][0] == 0:
                        #fullread mapped to genome but not in the same position with junction reads or shift of position on chromosome
                        #we remove these junction read as fullreads too
                    #    teFullReads[name] = 1
                elif jun_read_name[-2:] == '/2':
                    if not teJunctionReads.has_key(jun_read_name[:-2]) or not teJunctionReads[jun_read_name[:-2]].has_key(2):
                        continue
                    #read has label for read2, we use label for both paired or unpaired reads
                    #if teJunctionReads[jun_read_name[:-2]][2] == 1:
                    #    teFullReads[name] = 1
                    if teJunctionReads[jun_read_name[:-2]][2][3] == 1:
                        teFullReads[name] = 1
                    elif teJunctionReads[jun_read_name[:-2]][2][0] == chro and int(teJunctionReads[jun_read_name[:-2]][2][1]) == start:
                        #fullread overlap with junction read
                        if int(teJunctionReads[jun_read_name[:-2]][2][2]) >= end + extend:
                            #fullread extend more than 5 bp into insertion, suggesting no insertion
                            teFullReads[name] = 1
                    elif teJunctionReads[jun_read_name[:-2]][2][0] == chro and int(teJunctionReads[jun_read_name[:-2]][2][2]) == end:
                        if int(teJunctionReads[jun_read_name[:-2]][2][1]) <= start - extend:
                            teFullReads[name] = 1
                    #elif not teJunctionReads[jun_read_name[:-2]][2][0] == 0:
                    #    teFullReads[name] = 1
                elif read_order == 1 or read_order == 2:
                    if not teJunctionReads.has_key(jun_read_name) or not teJunctionReads[jun_read_name].has_key(read_order):
                        continue
                    #read do not have label, read is paired and read is read1
                    #if teJunctionReads[jun_read_name][1] == 1:
                        #fullread mapped to genome
                    #    teFullReads[name] = 1
                    if teJunctionReads[jun_read_name][read_order][3] == 1:
                        teFullReads[name] = 1
                    elif teJunctionReads[jun_read_name][read_order][0] == chro and int(teJunctionReads[jun_read_name][read_order][1]) == start:
                        #fullread overlap with junction read
                        if int(teJunctionReads[jun_read_name][read_order][2]) >= end + extend:
                            #fullread extend more than 5 bp into insertion, suggesting no insertion
                            teFullReads[name] = 1
                    elif teJunctionReads[jun_read_name][read_order][0] == chro and int(teJunctionReads[jun_read_name][read_order][2]) == end:
                        if int(teJunctionReads[jun_read_name][read_order][1]) <= start - extend:
                            teFullReads[name] = 1
                    #elif not teJunctionReads[jun_read_name][read_order][0] == 0:
                    #    teFullReads[name] = 1
                #elif read_order == 2:
                    #read do not have label, read is paired and read is read2
                #    if teJunctionReads[jun_read_name][2] == 1:
                        #fullread mapped to genome
                #        teFullReads[name] = 1
                else:
                    #read do not have label: read is unpaired, need to use unPaired_read_info
                    if teReadsInfo.has_key(name):
                        if teReadsInfo[name] == 1:
                            if not teJunctionReads.has_key(jun_read_name) and not teJunctionReads[jun_read_name].has_key(1):
                                continue
                            #read1
                            #if teJunctionReads[jun_read_name][1] == 1:
                            #    teFullReads[name] = 1
                            if teJunctionReads[jun_read_name][1][3] == 1:
                                teFullReads[name] = 1
                            elif teJunctionReads[jun_read_name][1][0] == chro and int(teJunctionReads[jun_read_name][1][1]) == start:
                                #fullread overlap with junction read
                                if int(teJunctionReads[jun_read_name][1][2]) >= end + extend:
                                    #fullread extend more than 5 bp into insertion, suggesting no insertion
                                    teFullReads[name] = 1
                            elif teJunctionReads[jun_read_name][1][0] == chro and int(teJunctionReads[jun_read_name][1][2]) == end:
                                if int(teJunctionReads[jun_read_name][1][1]) <= start - extend:
                                    teFullReads[name] = 1
                            #elif not teJunctionReads[jun_read_name][1][0] == 0:
                            #    teFullReads[name] = 1
                        elif teReadsInfo[name] == 2:
                            if not teJunctionReads.has_key(jun_read_name) and not teJunctionReads[jun_read_name].has_key(2):
                                continue
                            #read2
                            #if teJunctionReads[jun_read_name][2] == 1:
                            #    teFullReads[name] = 1
                            if teJunctionReads[jun_read_name][2][3] == 1:
                                teFullReads[name] = 1
                            elif teJunctionReads[jun_read_name][2][0] == chro and int(teJunctionReads[jun_read_name][2][1]) == start:
                                #fullread overlap with junction read
                                if int(teJunctionReads[jun_read_name][2][2]) >= end + extend:
                                    #fullread extend more than 5 bp into insertion, suggesting no insertion
                                    teFullReads[name] = 1
                            elif teJunctionReads[jun_read_name][2][0] == chro and int(teJunctionReads[jun_read_name][2][2]) == end:
                                if int(teJunctionReads[jun_read_name][2][1]) <= start - extend:
                                    teFullReads[name] = 1
                            #elif not teJunctionReads[jun_read_name][2][0] == 0:
                            #    teFullReads[name] = 1
                    else:
                        print 'teReadsInfo do not have key of %s' %(name)
                    
                #print '%s\t%s\t%s' %(name, jun_read_name, read_order)
                #if jun_read_name[-2:] == '/1':
                #    if teJunctionReads[jun_read_name[:-2]][1] == 1:
                #        teFullReads[name] = 1
                        #continue
                #elif jun_read_name[-2:] == '/2':
                #    if teJunctionReads[jun_read_name[:-2]][2] == 1:
                #        teFullReads[name] = 1
                        #continue
                #else:
                    #if teJunctionReads[jun_read_name][1] == 1 and teJunctionReads[jun_read_name][2] == 1:
                        #both reads are perfect matched
                        #in case the read is not paired because the pairs is in repeat
                    #    teFullReads[name] = 1
                        #continue
                    #elif teJunctionReads[jun_read_name][read_order] == 1:
                    #    teFullReads[name] = 1
                        #continue

            strand = ''
            cg_flag= 0
            # flag is 0 is read if read is unpaired and mapped to plus strand
            if int(flag) == 0:
                strand = '+'
            else:
                strand = '-' if record.is_reverse else '+'
            if r_cg1.search(cigar):
                #continue
                pass
            if r_cg2.search(cigar): 
                cg_flag = 1
            tags = convert_tag(tag)
            if verbose > 4: print '%s\t%s\t%s\t%s' %(name, start, MAPQ, tag)
            # filter low quality mapping reads: 
            # 1. paired-end reads at least have one reads unique mapped (MAPQ set to 0 for both reads if both are repeat, else should be > 0 at least one unique mapped)
            # 2. unpaired reads should unique mapped, no gap, mismatch <= 3 and no suboptimal alignment
            #print 'before: %s\t%s\t%s' %(name, count, bin_ins)
            #if record.is_proper_pair and (int(MAPQ) >= 29 or tags['XT'] == 'U'):
            #if record.is_proper_pair and int(MAPQ) > 0:
            if record.is_proper_pair:
                if verbose > 4: print 'is proper paired'
                #store low quality read in dict
                if int(MAPQ) < 29:
                    x1 = int(tags['X1']) if tags.has_key('X1') else 0
                    xm = int(tags['XM']) if tags.has_key('XM') else 0
                    xo = int(tags['XO']) if tags.has_key('XO') else 0
                    teLowQualityReads[name] = [int(MAPQ), xm, x1, xo]
                    
                if r.search(name): #junction reads, allowed 2 mismatch and 1 indels
                    if int(tags['XM']) <= int(mm_allow) and int(tags['XO']) <= 1:
                        bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, chro, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
                else:              #supporting reads, allowed upto 3 mismatch
                    if int(tags['XM']) <= int(mm_allow):
                        bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, chro, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
            #elif not record.is_paired:
            else:
                if verbose > 4: print 'is not paired or not proper paired'
                #store low quality read in dict, not properly paired is low quality
                if int(MAPQ) < 29 or record.is_paired:
                    if verbose > 4: print 'low quality reads'
                    x1 = int(tags['X1']) if tags.has_key('X1') else 0
                    xm = int(tags['XM']) if tags.has_key('XM') else 0
                    xo = int(tags['XO']) if tags.has_key('XO') else 0
                    teLowQualityReads[name] = [int(MAPQ), xm, x1, xo]
                #if tags['XT'] == 'U' and int(tags['XO']) == 0 and int(tags['XM']) <= 3 and int(tags['X1']) == 0:
                #if tags['XT'] == 'U' and int(tags['XO']) == 0 and (int(tags['XM']) <= 3 or int(tags['X1']) == 0):
                #if tags['XT'] == 'U' and int(tags['XO']) == 0 and int(tags['X1']) <= 3:
                if r.search(name): #junction reads, allowed only 1 mismatch
                    if tags['XT'] == 'U' and int(tags['XM']) <= int(mm_allow) and int(tags['X1']) <= 3:
                        if verbose > 4: print 'junction pass'
                        bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, chro, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
                else:
                    if tags['XT'] == 'U' and int(tags['XM']) <= int(mm_allow) and int(tags['X1']) <= 3:
                        if verbose > 4: print 'supporting pass'
                        bin_ins, count = align_process(bin_ins, read_repeat, record, r, r_tsd, count, seq, chro, start, end, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teReadClusters, teReadClusters_count, teReadClusters_depth, teSupportingReads)
            teReadClusters[count]['read_inf']['seq']['chr'] = chro
            #print 'after: %s\t%s\t%s' %(name, count, bin_ins)

    ###TSD not given we infer from read depth
    if r_tsd.search(TSD):
        TSD_from_read_depth(r, read_repeat, teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found, teLowQualityReads)        
        
def find_insertion_cluster_sam(sorted_align, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found):
    align_seg = pysam.AlignedSegment()
    r = re.compile(r':(start|end):(5|3)')
    r_tsd = re.compile(r'UNK|UKN|unknown', re.IGNORECASE)
    bin_ins        = [0]
    count          = 0
    TSD_len        = len(TSD)
    ##go through the reads, cluster reads and determine junction reads and supporting reads
    for record in sorted_align.keys():
        name, flag, ref, start, MAPQ, cigar, MRNM, MPOS, TLEN, seq, qual, tag = ['']*12
        try:
            name, flag, ref, start, MAPQ, cigar, MRNM, MPOS, TLEN, seq, qual, tag = re.split(r'\t', record, 11)
        except:
            try:
                name, flag, ref, start, MAPQ, cigar, MRNM, MPOS, TLEN, seq, qual = re.split(r'\t', record, 10)
            except:
                print 'sam format error, can not split into 11 or 12 field'
                exit()
        #tags = re.split(r'\t', tags)
        if int(MAPQ) < 40:
            continue
        end = len(seq) + int(start) - 1 # should not allowed for indel or softclip
        align_seg.flag = int(flag)
        strand = ''
        if int(flag) == 0:
            strand = '+'
        else:
            strand = '-' if align_seg.is_reverse else '+'
        #print name, TLEN, tags[0]
        range_allowance = 0
        padded_start    = bin_ins[0] - range_allowance
        padded_end      = bin_ins[-1] - range_allowance 
        m = r.search(name)
        if m:
            #insertions
            #print 'insertions: %s' %(name)
            if (int(start) >= padded_start and int(start) <= padded_end) or (int(end) >= padded_start and int(end) <= padded_end):
                bin_ins.extend([int(start), int(end)])
                bin_ins = sorted(bin_ins, key=int)
                if not r_tsd.search(TSD):
                    TSD_check(count, seq, start, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
                else:
                    calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
            else:
                #if start and end do not fall within last start and end
                #we now have a different insertion event
                count += 1
                if not r_tsd.search(TSD):
                    TSD_check(count, seq, start, name, TSD, strand, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)
                else:
                    calculate_cluster_depth(count, seq, start, name, strand, teReadClusters, teReadClusters_count, teReadClusters_depth)
                #initial insertion site boundary
                bin_ins = [int(start), int(end)]
        else:
            #supporting reads
            #print 'supportings: %s' %(name)
            pass
        #print name, flag, strand
    if r_tsd.search(TSD):
        TSD_from_read_depth(teReadClusters, teReadClusters_count, teReadClusters_depth, teInsertions, teInsertions_reads, existingTE_inf, existingTE_found)

##default sam, not deal with paired end mapping and supporting reads 
def remove_redundant_sam(infile, target):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'@'):
                unit = re.split(r'\t',line)
                if target != 'ALL' and unit[2] != target:
                    continue
                else:
                    start = unit[3]
                    if not data.has_key(line):
                        data[line] = start
    data = OrderedDict(sorted(data.items(), key=lambda x: x[1]))
    return data

def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def read_junction_reads_align(align_file_f, read_repeat, teJunctionReads):
    target = 'ALL'
    ref  = None if target == 'ALL' else target
    fsam = pysam.AlignmentFile(align_file_f, 'rb')
    rnames = fsam.references
    for record in fsam.fetch(reference=ref, until_eof = True):
        if record.query_sequence is None:
            continue
        read_order = 0
        if record.is_read1:
            read_order = 1
        elif record.is_read2:
            read_order = 2 
        if not record.is_unmapped:
            name   = record.query_name
            flag   = record.flag
            MAPQ   = record.mapping_quality
            cigar  = record.cigarstring
            seq    = record.query_sequence
            tag    = record.tags if record.tags else []
            chro   = rnames[record.reference_id]
            start  = int(record.reference_start) + 1
            end    = int(record.reference_end) + 1
            match  = 0
            tags   = convert_tag(tag)
            if verbose > 4: print 'read junction align: %s\t%s\t%s\t%s' %(name, read_order, seq, tag)
            for (key, length) in record.cigartuples:
                #print key, length
                if int(key) == 0:
                    match += length
                elif int(key) == 1:
                    #treat insertion in reads as match, which allow indel in alignment
                    match += length
            #if int(MAPQ) >= 29 and match >= len(seq) - 10:
            #    print '%s' %(name)
            intact_flag = 0
            if match >= len(seq) - 10:
                #teJunctionReads[name][read_order] = 1
                intact_flag = 1
          
            ##quality control of fullreads alignment
            #xm = int(tags['XM']) if tags.has_key('XM') else 0
            #if xm <= 3:
            read_name0 = name
            read_name1 = '%s/1' %(name)
            read_name2 = '%s/2' %(name)
            if read_repeat.has_key(read_name0) or read_repeat.has_key(read_name1) or read_repeat.has_key(read_name2):
                teJunctionReads[name][read_order] = [chro, start, end, intact_flag]
                if verbose > 4: print 'TE_junction_reads, chromosome specific: %s\t%s\t%s\t%s' %(name, chro, start, end)
            #else:
            #    teJunctionReads[name][read_order] = [0, 0, 0, 0]
        else:
            name   = record.query_name
            read_name0 = name
            read_name1 = '%s/1' %(name)
            read_name2 = '%s/2' %(name)
            if read_repeat.has_key(read_name0) or read_repeat.has_key(read_name1) or read_repeat.has_key(read_name2):
                teJunctionReads[name][read_order] = [0, 0, 0, 0]

def read_unpaired_read_info(unpaired_read_info, teReadsInfo):
    with open (unpaired_read_info, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                teReadsInfo[unit[0]] = int(unit[1])

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
    if not len(sys.argv) == 13:
        usage()
        exit(2)

    global verbose
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
    verbose              = int(sys.argv[12])
    #relax_reference      = sys.argv[11]## relax mode for existing TE: 1 or 0
    #relax_align          = sys.argv[12]## relax mode for insertion: 1 or 0
    bowtie_sam           = 1           ## change to shift or remove in V2
    existingTE_inf       = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str)))
    existingTE_found     = defaultdict(lambda : defaultdict(lambda : int))
    bwa                  = 0
    teInsertions         = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int()))))
    teInsertions_reads   = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list()))))
    teReadClusters       = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str()))))
    teReadClusters_count = defaultdict(lambda : defaultdict(lambda : int()))
    teReadClusters_depth = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int()))))
    teSupportingReads    = defaultdict(lambda : list())
    teLowQualityReads    = defaultdict(lambda : list())
    teFullReads          = defaultdict(lambda : int())
    teJunctionReads      = defaultdict(lambda : defaultdict(lambda : list()))
    teReadsInfo          = defaultdict(lambda : int())

    top_dir = re.split(r'/', os.path.dirname(os.path.abspath(align_file)))[:-1]
    #read existing TE from file
    r_te = re.compile(r'repeatmasker|rm|\.out', re.IGNORECASE)
    if os.path.isfile(existing_TE) and os.path.getsize(existing_TE) > 0:
        if r_te.search(existing_TE):
            existingTE_RM_ALL('/'.join(top_dir), existing_TE, existingTE_inf, usr_target)
        #else:
        #    existingTE('/'.join(top_dir), existing_TE, existingTE_inf, existingTE_found)
    else:
        print 'Existing TE file does not exists or zero size'

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
    read_repeat = read_repeat_name(read_repeat_files)

    #read and store full length reads alignment on genome
    #*.repeat.fullreads.bwa.sorted.bam
    #*.repeat.bwa.sorted.bam
    align_file_f = re.sub(r'.bwa.sorted.bam', r'.fullreads.bwa.sorted.bam', align_file)
    print 'fullread bam: %s' %(align_file_f)
    read_junction_reads_align(align_file_f, read_repeat, teJunctionReads)
    print 'junction fullread number: %s' %(len(teJunctionReads.keys()))   
 
    #read and store unpaired junction read information: read1 or read2
    unpaired_read_info = '%s/flanking_seq/%s.flankingReads.unPaired.info' %('/'.join(top_dir), usr_target)
    read_unpaired_read_info(unpaired_read_info, teReadsInfo) 

    ##cluster reads around insertions
    find_insertion_cluster_bam(align_file, read_repeat, usr_target, TSD, teInsertions, teInsertions_reads, teReadClusters, teReadClusters_count, teReadClusters_depth, existingTE_inf, existingTE_found, teSupportingReads, teLowQualityReads, teJunctionReads, teFullReads, teReadsInfo, mm_allow)


    ##output insertions
    write_output(top_dir, result, read_repeat, usr_target, exper, TE, required_reads, required_left_reads, required_right_reads, teInsertions, teInsertions_reads, teSupportingReads, existingTE_inf, teReadClusters, bedtools, lib_size, teLowQualityReads, teFullReads)

 
if __name__ == '__main__':
    main()

