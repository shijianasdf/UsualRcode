# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

# 基础设置
export xteam=/pub5/xiaoyun/Bin

export PANDOC_HOME=/pub5/xiaoyun/Software/pandoc-2.1.1
export GCC_HOME=/IJob/J21/gcc-7.1.0
export texlive=/usr/local/texlive/2017/bin/x86_64-linux
#special for GISTIC2.0
export mcr_root=/pub5/xiaoyun/BioSoftware/MATLAB/MATLAB_Compiler_Runtime
export XAPPLRESDIR=$mcr_root/v714/X11/app-defaults

# JAVA
export JAVA_HOME=/usr/lib/jvm/java-1.7.0-openjdk.x86_64
export CLASSPATH=$JAVA_HOME/lib:$JAVA_HOME/jre/lib:$CLASSPATH

# 比对及装配
export Samtools=/pub5/xiaoyun/BioSoftware/Samtools/bin
export bowtie=/pub5/xiaoyun/BioSoftware/bowtie
export bowtie_index=/pub5/xiaoyun/BioSoftware/bowtie/indexes/hg19
export bowtie_threads=24
export bowtie2=/pub5/xiaoyun/BioSoftware/bowtie2
export tophat2=/pub5/xiaoyun/BioSoftware/tophat2
export cufflinks=/pub5/xiaoyun/BioSoftware/Cufflinks

# sra2fastq
export sratoolkit=/pub5/xiaoyun/BioSoftware/sratoolkit.2.7/bin

# 质量控制
export FastQC=/pub5/xiaoyun/BioSoftware/FastQC

# bedtools
export bedtools=/pub5/xiaoyun/BioSoftware/bedtools2

# 下载器
export ASPERA=/pub5/xiaoyun/Software/aspera/connect

# perl
export PERL_HOME=/pub5/xiaoyun/BioSoftware/perl-5.22.0

# python环境变量
# export PYTHON_HOME=/pub5/xiaoyun/BioSoftware/Python-2.7.12
export PYTHON_HOME=/pub5/xiaoyun/BioSoftware/Python-3.6.5

# motif 分析
export MEMEHOME=/pub5/xiaoyun/BioSoftware/meme-4.11.2
export HOMER_HOME=/pub5/xiaoyun/BioSoftware/HOMER
export WEBLOGO=/pub5/xiaoyun/BioSoftware/weblogo
export MPIHOME=/pub5/xiaoyun/BioSoftware/openmpi-2.0.1
export ENCODEMOTIF_HOME=/pub5/xiaoyun/BioSoftware/motifAnalysis/encodeMotifs

# UCSC 各种可执行文件
export UCSC_HOME=/pub5/xiaoyun/BioSoftware/ucscSuite

# VPN代理
export PROXYCHAINS_CONF_FILE=/pub5/xiaoyun/BioSoftware/proxychains/etc/proxychains.conf
export PROXYCHAINS_HOME=/pub5/xiaoyun/BioSoftware/proxychains

# R 
export pcre_PATH=/pub5/xiaoyun/Software/pcre-8.42
# export XZ_HOME=/pub5/xiaoyun/Software/R-replay/xz-5.2.4
export XZ_HOME=/pub5/xiaoyun/Software/xz-5.2.3
# export XZ_HOME=/pub5/xiaoyun/Software/R-replay/xz-5.0.8



# 组合LD_LIBRARY_PATH：为了保证新添加的能够覆盖系统默认，需要把新的地址放在前面，系统默认放在最后
export LD_LIBRARY_PATH=$GCC_HOME/lib64:$mcr_root/v714/runtime/glnxa64:$mcr_root/v714/sys/os/glnxa64:$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64:$XZ_HOME/lib:/usr/lib:/usr/local/lib:/usr/lib64:/usr/local/gmp-6.1.0/lib:/usr/local/mpc-1.0.3/lib:/usr/local/mpfr-3.1.4/lib:$pcre_PATH/lib::/usr/lib:/usr/local/lib:/usr/lib64$LD_LIBRARY_PATH


# 组合PATH
export PATH=$xteam:$PANDOC_HOME/bin:$GCC_HOME/bin:$GCC_HOME:$texlive:$XZ_HOME/bin:$JAVA_HOME:$JAVA_HOME/jre:$Samtools:$bowtie:$bowtie2:$tophat2:$cufflinks:$sratoolkit:$FastQC:$bedtools/bin:$ASPERA/bin:$PERL_HOME/bin:$PYTHON_HOME:$PYTHON_HOME/bin:$MEMEHOME/bin:$HOMER_HOME/bin:$WEBLOGO:$MPIHOME/bin:$ENCODEMOTIF_HOME/encode-motifs-v1.3/bin:$UCSC_HOME:$PROXYCHAINS_HOME/bin:$pcre_PATH/bin:$PATH

#gd
#export GDHOME=/pub5/xiaoyun/Software/gd-2.0.28
#PATH=$GDHOME/bin:$PATH



# 组长特异的
PERL_MB_OPT="--install_base \"/home/xiaoyun/perl5\""; export PERL_MB_OPT;
PERL_MM_OPT="INSTALL_BASE=/home/xiaoyun/perl5"; export PERL_MM_OPT;

