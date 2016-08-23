echo "backup CONFIG file"
backup_time=`date | sed 's/ \|:/_/g'`
mv CONFIG CONFIG\.$backup_time
echo "#tools" > CONFIG
config=`pwd`/CONFIG
echo "" > install.log
log=`pwd`/install.log
pythonlib=`pwd`/pythonlib
bin=`pwd`/bin

if [ ! -d './pythonlib' ]; then
   mkdir pythonlib
fi

if [ ! -d './bin' ]; then
   mkdir bin
fi

echo "installing external tools into External_tools"
if [ ! -d './External_tools' ]; then
    mkdir External_tools
fi
cd External_tools

#echo "installing bwa-0.6.2"
if [ ! -f './bwa-0.6.2/bwa' ]; then
    echo "installing bwa-0.6.2"
    wget --output-document bwa-0.6.2.tar.bz2 'http://downloads.sourceforge.net/project/bio-bwa/bwa-0.6.2.tar.bz2?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbio-bwa%2Ffiles%2F&ts=1471890846&use_mirror=superb-sea2' >> $log 2>&1
    tar -jxf bwa-0.6.2.tar.bz2 >> $log 2>&1
    cd bwa-0.6.2
    make >> $log 2>&1
    if test -e bwa; then
        echo "bwa installed"
        #echo "bwa=`pwd`/bwa" >> $config
    else
        echo "bwa installation failed. Try to install manually."
    fi
    cd ..
fi

#echo "installing seqtk released 1.2"
if [ ! -f './seqtk/seqtk' ]; then
    echo "installing seqtk released 1.2"
    git clone https://github.com/lh3/seqtk.git >> $log 2>&1
    cd seqtk
    make >> $log 2>&1
    if test -e seqtk; then
        echo "seqtk installed"
        #echo "seqtk=`pwd`/seqtk" >> $config
    else
        echo "seqtk installation failed. Try to install manually."
    fi
    cd ..
fi

#echo "installing bowtie2.2.9"
if [ ! -f './bowtie2-2.2.9/bowtie2' ]; then
    echo "installing bowtie2.2.9"
    wget --output-document bowtie2-2.2.9-linux-x86_64.zip 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie2%2F2.2.9%2F&ts=1471884842&use_mirror=pilotfiber' >> $log 2>&1
    unzip bowtie2-2.2.9-linux-x86_64.zip >> $log 2>&1
    cd bowtie2-2.2.9
    if test -e bowtie2; then
        echo "bowtie2 installed"
        #echo "bowtie2=`pwd`/bowtie2" >> $config
    else
        echo "bowtie2 installation failed. Try to install manually."
    fi
    cd ..
fi

#echo "installing blat"
if [ ! -f './blat/blat' ]; then
    echo "installing blat"
    mkdir blat
    cd blat
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat >> $log 2>&1
    chmod 755 blat
    if test -e blat; then
        echo "blat installed"
        #echo "blat=`pwd`/blat" >> $config
    else
        echo "blat installation failed. Try to install manually."
    fi
    cd .. 
fi

#echo "installing pysam"
if [ ! -d $pythonlib/lib64/python2.7/site-packages/pysam-0.9.1.4-py2.7-linux-x86_64.egg ]; then
    echo "installing pysam"
    git clone https://github.com/pysam-developers/pysam.git >> $log 2>&1
    cd pysam
    export PYTHONPATH=$pythonlib/lib64/python2.7/site-packages:$PYTHONPATH
    mkdir -p $pythonlib/lib64/python2.7/site-packages
    python setup.py install --prefix $pythonlib >> $log 2>&1
    cd ../..
    (PYTHONPATH=$pythonlib/lib64/python2.7/site-packages python -c "import pysam;print pysam.__version__" && echo "pysam installed" ) || echo "pysam installation failed. Try to install manually."
    cd External_tools
fi

#echo "installing bedtools"
if [ ! -f './bedtools2/bin/bedtools' ]; then
    echo "installing bedtools"
    wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz >> $log 2>&1
    tar -zxvf bedtools-2.25.0.tar.gz >> $log 2>&1
    cd bedtools2
    make >> $log 2>&1
    if test -e ./bin/bedtools; then
        echo "bedtools installed"
        #echo "bedtools=`pwd`/bin/bedtools" >> $config
    else
        echo "bedtools installation failed. Try to install manually."
    fi
    cd ..
fi

#echo "installing samtools"
if [ ! -f './samtools-1.3.1/samtools' ]; then
    echo "installing samtools"
    wget --output-document samtools-1.3.1.tar.bz2 'http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Fsamtools%2F1.3.1%2F&ts=1471888881&use_mirror=heanet' >> $log 2>&1
    tar -jxf samtools-1.3.1.tar.bz2 >> $log 2>&1
    cd samtools-1.3.1
    ./configure --enable-plugins --enable-libcurl --prefix=`pwd` >> $log 2>&1
    make all all-htslib >> $log 2>&1
    make install install-htslib >> $log 2>&1
    if test -e samtools; then
        echo "samtools installed"
        #echo "samtools=`pwd`/samtools" >> $config
    else
        echo "samtools installation failed. Try to install manually."
    fi
    cd ..
fi
#back to RelocaTE2 home directory
cd ..

test_exe ()
{
    #echo $1
    #echo $2
    if test -e $1; then
        echo "$3 installed"
        if [ ! -e $2 ]; then
            ln -s $1 $2
        fi
        echo "$3=$2" >> $4
    else
        echo "$3 installation failed. Try to install manually."
    fi
}

echo ""
echo "Testing installation and making CONFIG file"

test_exe `pwd`/External_tools/bwa-0.6.2/bwa $bin/bwa bwa $config
test_exe `pwd`/External_tools/seqtk/seqtk $bin/seqtk seqtk $config
test_exe `pwd`/External_tools/bowtie2-2.2.9/bowtie2 $bin/bowtie2 bowtie2 $config
test_exe `pwd`/External_tools/bowtie2-2.2.9/bowtie2-build $bin/bowtie2-build bowtie2_build $config
test_exe `pwd`/External_tools/blat/blat $bin/blat blat $config
test_exe `pwd`/External_tools/bedtools2/bin/bedtools $bin/bedtools bedtools $config
test_exe `pwd`/External_tools/samtools-1.3.1/samtools $bin/samtools samtools $config

echo " "
echo "New CONFIG file is:"
cat $config
echo " "

echo "Done"

