####bedtools根据bed文件找到fasta文件,bed文件一般不要表头####
/pub5/xiaoyun/BioSoftware/bedtools2/bin/bedtools getfasta -fi /pub5/xiaoyun/BioSoftware/bowtie/indexes/hg19.fa -bed region.bed -fo region.fa
/pub5/xiaoyun/BioSoftware/bedtools2/bin/bedtools getfasta -fi /pub5/xiaoyun/BioSoftware/bowtie/indexes/hg19.fa -bed region.bed -fo region.fa
###bedtools基于染色体进行排序,bed文件一般不要表头,进行merge###
sortBed -i super-enhancer.bed > super-enhancer-sort.bed #sort后的顺序不是按chr1,chr2....这个顺序，而是按chr1,chr10,chr11.....chr19,chr2,chr20..chr22,chr3
bedtools merge -i super-enhancer-sort.bed > super-enhancer-sort-merge.bed

####bedtools求两个区间集合的overlap,首先bed文件是没有头的,一般都是先sortBed####
/pub5/xiaoyun/BioSoftware/bedtools2/bin/sortBed  -i /home/shijian2015/geneHancer.region.bed > /home/shijian2015/geneHancer.region.sort.bed 
/pub5/xiaoyun/BioSoftware/bedtools2/bin/sortBed  -i /home/shijian2015/TFBS.region.bed > /home/shijian2015/TFBS.region.sort.bed 
/pub5/xiaoyun/BioSoftware/bedtools2/bin/intersectBed -a /home/shijian2015/geneHancer.region.sort.bed -b /home/shijian2015/TFBS.region.sort.bed -wa -wb  > /home/shijian2015/result.bed
ClosestBed  #找最近的基因

#####bedtools merge,先sort再merge#####
/pub5/xiaoyun/BioSoftware/bedtools2/bin/sortBed  -i /home/shijian2015/geneHancer.region.bed > /home/shijian2015/geneHancer.region.sort.bed 
/pub5/xiaoyun/BioSoftware/bedtools2/bin/mergeBed -i /pub6/temp/shijian/enhancer.TFBS.overlap/geneHancer.region.sort.bed   > /pub6/temp/shijian/enhancer.TFBS.overlap/geneHancer.region.sort.merge.bed

####transfac2meme根据transfac motif文件得到meme需要的motif文件格式####
transfac2meme enrichLncBsEnmotif.dat > enrichLncBsEnmotif.meme

################################################################FIMO#############################################################
#FIMO:返回TF motif在一个序列或一堆序列上匹配的位置
#FIMO converts each input motif into a log-odds PSSM and uses each PSSM to independently scan each input sequence. It reports all positions in each sequence that match a motif with a statistically significant log-odds score. You can control the match p-value that is considered significant, and whether or not FIMO reports matches on both strands when the sequence alphabet is complementable (e.g., DNA or RNA).
#FIMO网页版输入meme motif数据格式(pwm转换而来)，以及fasta格式数据
#Usage:fimo [options] <motifs> <database>  
#motifs:A file containing MEME formatted motifs. MEME格式的motif,可以用MEME suite自带的脚本程序去把其他的motif格式转换成MEME格式的motif
#database:A file containing a collection of sequences in FASTA format. FASTA格式的序列集合
#options:各种参数的设定 -o 指定输出目录
#meme motif格式:
##################################################################################################################################
/pub5/xiaoyun/BioSoftware/meme-4.11.2/bin/fimo -o /pub6/temp/shijian/enhancer.TFBS.overlap/result /pub6/temp/shijian/enhancer.TFBS.overlap/MA0014.3.meme /pub6/temp/shijian/enhancer.TFBS.overlap/top10_geneHancer.region.fa

################################################################MEME#############################################################
#MEME:discovers novel, ungapped motifs (recurring, fixed-length patterns) in your sequences
#MEME takes as input a group of sequences and outputs as many motifs as requested
#Usage: meme <sequence file> [options]
#sequence file:The name of a file containing FASTA sequences or the word stdin to indicate that the sequences should be read from standard input.
#options:各种参数的设定 -o 指定输出目录 -dna 指定DNA链
##################################################################################################################################
export MEME_HOME=/pub5/xiaoyun/BioSoftware/meme-4.11.2/bin
export PATH=$PATH:$MEME_HOME/bin
meme -dna /pub6/temp/shijian/enhancer.TFBS.overlap/top10_geneHancer.region.fa -o /pub6/temp/shijian/enhancer.TFBS.overlap/meme.result


################HOMER进行motif富集分析(HOMER既可以进行denove motif discovery analysis也可以进行motif富集分析)######################
# 设置HOMER环境变量
export HOMER_HOME=/pub5/xiaoyun/BioSoftware/HOMER
export WEBLOGO_HOME=/pub5/xiaoyun/BioSoftware/weblogo
export PATH=$PATH:$HOMER_HOME/bin:$WEBLOGO_HOME
findMotifsGenome.pl /IJob/J26/lunwen/3.characterize.cancer.associated.enhacner/SNPenrichment/LUSC/LUSC_specific_enhancers_region_hg19_homer.bed hg19 /IJob/J26/lunwen/3.characterize.cancer.associated.enhacner/motif/LUSC_homer.result -size 200 -mask -basic -fdr 0.05
findMotifsGenome.pl /IJob/J26/lunwen/3.characterize.cancer.associated.enhacner/SNPenrichment/breast/breast_specific_enhancers_region_hg19_homer.bed hg19 /IJob/J26/lunwen/3.characterize.cancer.associated.enhacner/motif/BRCA_homer.result -size 200 -mask -basic -fdr 0.05

#############liftover:基因组坐标转换################
cd /pub5/xiaoyun/Jobs/J23/TAD_Disruption/Software
./liftOver /pub5/xiaoyun/Jobs/J23/TAD_Disruption/Data/GSE35156/hESC_domain/combined/total.combined.domain  /pub6/temp/shijian/hg19ToHg38.over.chain.gz  /pub5/xiaoyun/Jobs/J23/TAD_Disruption/Data/GSE35156/hESC_domain/combined/total.combined.domain_hg19 unlifted.bed

##################aspera下载数据命令行#################
/pub5/xiaoyun/Software/aspera/connect/bin/ascp -i /pub5/xiaoyun/Software/aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l300m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR713/SRR7138443/SRR7138443.sra /IData/CancerOMICS/BrainTumour\(SRP145073\)/ 

##################rsync同步数据####################
rsync -aP rsync://ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR401/SRR4019253/SRR4019253.sra FU97/GSM2267399/1.SRA/





