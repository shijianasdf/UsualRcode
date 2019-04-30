# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs
#设置环境变量
PATH=$PATH:$HOME/bin

export PATH

export xteam=/pub5/xiaoyun

#查看变量
echo $xteam

# 测序数据相关的基础软件
# 比对及装配
export bowtie=/pub5/xiaoyun/BioSoftware/bowtie
export bowtie2=/pub5/xiaoyun/BioSoftware/bowtie2
export tophat2=/pub5/xiaoyun/BioSoftware/tophat2
export cufflinks=/pub5/xiaoyun/BioSoftware/Cufflinks
# samtools
export Samtools=/pub5/xiaoyun/BioSoftware/Samtools/bin/
# sra2fastq
#export sratoolkit=/pub5/xiaoyun/BioSoftware/sratoolkit/bin/
export sratoolkit=/pub5/xiaoyun/BioSoftware/sratoolkit.2.7/bin/
# 质量控制
export FastQC=/pub5/xiaoyun/BioSoftware/FastQC
# bedtools
export bedtools=/pub5/xiaoyun/BioSoftware/bedtools2
# 下载器
export ASPERA=/pub5/xiaoyun/Jobs/temp/aspera/connect
# perl
export PERL_HOME=/pub5/xiaoyun/BioSoftware/perl-5.22.0
# python环境变量
#export PYTHON27HOME=/pub5/xiaoyun/BioSoftware/Python-2.7.10
export PYTHONHOME=/pub5/xiaoyun/BioSoftware/Python-2.7.12
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/python2.6
# motif 分析套件
export MEMEHOME=/pub5/xiaoyun/BioSoftware/meme-4.11.2
#/pub5/xiaoyun/BioSoftware/motifAnalysis/meme.3.5.7
export TRAWLERHOME=/pub5/xiaoyun/BioSoftware/motifAnalysis/Trawler
export WEEDERHOME=/pub5/xiaoyun/BioSoftware/motifAnalysis/Weeder1.4.2
export ALIGNACEHOME=/pub5/xiaoyun/BioSoftware/motifAnalysis/alignace
export MDSCANHOME=/pub5/xiaoyun/BioSoftware/motifAnalysis/MDscan.2004
export ENCODEMOTIF=/pub5/xiaoyun/BioSoftware/motifAnalysis/encodeMotifs
export HOMER_HOME=/pub5/xiaoyun/BioSoftware/HOMER
export WEBLOGO_HOME=/pub5/xiaoyun/BioSoftware/weblogo
# UCSC 各种可执行文件
export ucscSuite=/pub5/xiaoyun/BioSoftware/ucscSuite

# 此部分专用于 chirpseq pipeline
export chirpPipe=/pub5/xiaoyun/BioSoftware/ChIRP.pipeline
# chirpseq pipeline 目录
export scripts_directory=/pub5/xiaoyun/BioSoftware/ChIRP.pipeline
# 染色体大小
export chromosome_size_file=/pub5/xiaoyun/BioSoftware/bedtools2/genomes/human.hg19.genome
# Bowtie 相关路径
export bowtie_bin=/pub5/xiaoyun/BioSoftware/bowtie/bowtie
export bowtie_index=/pub5/xiaoyun/BioSoftware/bowtie/indexes/hg19
export bowtie_threads=24
# MACS 路径
export macs_bin=/pub5/xiaoyun/BioSoftware/MACS14/bin/macs
# UCSC 程序路径
export wigToBigWig=/pub5/xiaoyun/BioSoftware/ucscSuite/wigToBigWig
export bigWigToBedGraph=/pub5/xiaoyun/BioSoftware/ucscSuite/bigWigToBedGraph
export bedGraphToBigWig=/pub5/xiaoyun/BioSoftware/ucscSuite/bedGraphToBigWig
export PROXYCHAINS_CONF_FILE=/pub5/xiaoyun/BioSoftware/proxychains/etc/proxychains.conf
export PROXYCHAINS=/pub5/xiaoyun/BioSoftware/proxychains

# 修改环境变量PATH
export PATH=$PROXYCHAINS/bin:$PERL_HOME/bin:$PYTHONHOME:$PYTHONHOME/bin:$bowtie:$bowtie2:$FastQC:$GISTIC2:$Samtools:$sratoolkit:$tophat2:$ASPERA/bin:$bedtools/bin:$cufflinks:$xteam/Bin:$MEME/bin:$ENCODEMOTIF/encode-motifs-v1.3:$ENCODEMOTIF/encode-motifs-v1.3/bin:$gettext/bin:$MEMEHOME/bin:$TRAWLERHOME/bin:$WEEDERHOME:$ALIGNACEHOME:$MDSCANHOME:$ucscSuite:$chirpPipe:$HOMER_HOME/bin:$WEBLOGO_HOME:$PATH
