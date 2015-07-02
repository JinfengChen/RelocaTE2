#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python Python.py --input ALL.all_nonref_insert.gff --refte ~/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa.RepeatMasker.out.bed

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

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

#Chr3	not.give	RelocaTE_i	283493	283504	.	-	.	ID=repeat_Chr3_283493_283504;TSD=ATGCCATCAAGG;Note=Non-reference,
#not found in reference;Right_junction_reads:4;Left_junction_reads:1;Right_support_reads:4;Left_support_reads:5;
#Chr3	281479	284272	TE110112	+
def Overlap_TE_boundary(prefix, refte, distance):
    data = defaultdict(str)
    final_gff   = '%s.gff' %(prefix)
    raw_gff     = '%s.raw.gff' %(prefix)
    all_gff     = '%s.all.gff' %(prefix)
    high_gff    = '%s.high_conf.gff' %(prefix)
    clean_gff   = '%s.clean.gff' %(prefix)
    infile = '%s.overlap' %(prefix)
    outfile= '%s.remove.gff' %(prefix)
    os.system('/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools window -w %s -a %s -b %s > %s' %(int(distance) + 10, final_gff, refte, infile))
    if not os.path.isfile(infile) or not os.path.getsize(infile) > 0:
        return 1
    print 'filter existing TE: %s bp' %(distance) 
    ofile  = open(outfile, 'w') 
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                temp = defaultdict(str)
                attrs = re.split(r';', unit[8])
                print line
                for attr in attrs:
                    if not attr == '':
                        attr = re.sub(r':', '=', attr)
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                if int(temp['Right_junction_reads']) == 0 or int(temp['Left_junction_reads']) == 0:
                    #support by one junction
                    #within 10 bp interval of intact TE boundary
                    #print >> ofile, '\t'.join(unit[:9])
                    if int(unit[3]) >= int(unit[10]) - int(distance) and int(unit[3]) <= int(unit[10]) + int(distance):
                        print >> ofile, '\t'.join(unit[:9])
                    elif int(unit[3]) >= int(unit[11]) - int(distance) and int(unit[3]) <= int(unit[11]) + int(distance):
                        print >> ofile, '\t'.join(unit[:9])
                    elif int(unit[4]) >= int(unit[10]) - int(distance) and int(unit[4]) <= int(unit[10]) + int(distance):
                        print >> ofile, '\t'.join(unit[:9])
                    elif int(unit[4]) >= int(unit[11]) - int(distance) and int(unit[4]) <= int(unit[11]) + int(distance):
                        print >> ofile, '\t'.join(unit[:9])
    ofile.close()
    if not os.path.isfile(outfile) or not os.path.getsize(outfile) > 0:
        #nothing to remove
        print 'nothing to remove'
        os.system('mv %s %s' %(final_gff, raw_gff))
        os.system('cp %s %s' %(raw_gff, all_gff))
        os.system('grep -v \"singleton\|insufficient_data\|supporting_reads\" %s > %s' %(raw_gff, final_gff))
        os.system('grep -v -e \"Right_junction_reads=1;Left_junction_reads=0\" -e \"Right_junction_reads=0;Left_junction_reads=1\" %s > %s' %(final_gff, high_gff)) 
        os.system('rm %s.overlap %s.remove.gff' %(prefix, prefix))
    else:
        print 'remove by bedtool'
        os.system('/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools intersect -v -a %s -b %s > %s' %(final_gff, outfile, clean_gff))
        os.system('mv %s %s' %(final_gff, raw_gff))
        os.system('cp %s %s' %(clean_gff, all_gff))
        os.system('grep -v \"singleton\|insufficient_data\|supporting_reads\" %s > %s' %(clean_gff, final_gff))
        os.system('grep -v -e \"Right_junction_reads=1;Left_junction_reads=0\" -e \"Right_junction_reads=0;Left_junction_reads=1\" %s > %s' %(final_gff, high_gff))
        os.system('rm %s.overlap %s.remove.gff %s.clean.gff' %(prefix, prefix, prefix))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--refte')
    parser.add_argument('-d', '--distance', default='3', type=int)
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.refte) > 0
    except:
        usage()
        sys.exit(2)
    

    #if not os.path.isfile(args.input) or not os.path.getsize(args.input) > 0:
    #txt2gff('%s.txt' %(os.path.splitext(args.input)[0]), args.input, 'non_reference')
    #os.system('bedtools window -w 10 -a %s -b %s > %s.overlap' %(args.input, args.refte, os.path.splitext(args.input)[0]))
    Overlap_TE_boundary(os.path.splitext(args.input)[0], args.refte, args.distance)
    

if __name__ == '__main__':
    main()

