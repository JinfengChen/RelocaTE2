#utility.py
#functions used in several scripts
#jinfeng.chen@ucr.edu



##############################dirs and files##################
##create new directory
def createdir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

##write content to new file
def writefile(outfile, lines):
    ofile = open(outfile, 'w')
    print >> ofile, lines
    ofile.close()


###############################multiprocess#####################
##run function with parameters using multiprocess of #cpu
def mp_pool_function(function, parameters, cpu):
    pool = mp.Pool(int(cpu))
    imap_it = pool.map(function, tuple(parameters))
    collect_list = []
    for x in imap_it:
        print 'status: %s' %(x)
        collect_list.append(x)
    return collect_list

##run command line by os.system
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



##########################fasta##############################
##get fasta id
def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

##complement sequence
def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    for i in range(len(bases)):
        bases[i] = complement[bases[i]] if complement.has_key(bases[i]) else bases[i]
    return ''.join(bases)

##reverse_complement sequence
def reverse_complement(seq):
    return complement(seq[::-1])



#########################bam################################
##convert tags collumn in sam format into dictionary
def convert_tag(tag):
    tags = {}
    for t in tag:
        tags[t[0]] = t[1]
    return tags



#########################relocate###########################
##parse regex.txt
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


