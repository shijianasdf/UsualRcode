
docker pull ubuntu
docker pull bcgsc/orca
docker pull austins2/cellranger6.0.2
docker run -it centos/httpd-24-centos7 /bin/bash
docker run -it conda/miniconda3-centos7 /bin/bash
docker run -it austins2/cellranger6.0.2 /bin/bash
exit
docker images
Docker search http
Docker search condon
Docker search anaconda
docker run -it   unlhcc/cellranger  /bin/bash
docker run -it -v /Users/shijian/mydata:/mnt/mydata  austins2/cellranger6.0.2  /bin/bash
docker run -it -v /Users/shijian/mydata:/mnt/mydata  unlhcc/cellranger  /bin/bash
docker run -it -v /Volumes/A/tianjin_mouse_single_cell:/mnt/tianjin_mouse_single_cell  austins2/cellranger6.0.2  /bin/bash

#cellRanger count
cd  /mnt/mydata/P2_656098
cellranger count --id=mRNA  \
--fastqs=/mnt/mydata/P2_656098/fastq/X101SC21081645-Z01-J004 \
--sample=P2_656098-1,P2_656098-2,P2_656098-3,P2_656098-4,P2_656098-5,P2_656098-6,P2_656098-7,P2_656098-8  \
--transcriptome=/mnt/mydata/refdata-gex-GRCh38-2020-A  

#tcr
cd /mnt/mydata/P2_656098
cellranger vdj --id=TCR \
               --fastqs=/mnt/mydata/P2_656098/fastq/X101SC21081645-Z01-J005 \
               --sample=P2_656098\\\
               --reference=/mnt/mydata/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
               --chain TR

#bcr
cd /mnt/mydata/P2_656098
cellranger vdj --id=BCR \
               --fastqs=/mnt/mydata/P2_656098/fastq/X101SC21081645-Z01-J006 \
               --sample=P2_656098 \\\
               --reference=/mnt/mydata/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
               --chain IG

#cellRanger count 
cd  /mnt/mydata/kidney/K1-10X5
cellranger count --id=mRNA   \
--fastqs=/mnt/mydata/kidney/LBFC20211146/210902_A00838_0544_BHKFJCDSX2/K1  \
--sample=K1-10X5   \
--transcriptome=/mnt/mydata/refdata-gex-GRCh38-2020-A  

cd  /mnt/mydata/kidney/K2-10X5
cellranger count --id=mRNA   \
--fastqs=/mnt/mydata/kidney/LBFC20211146/210902_A00838_0544_BHKFJCDSX2/K2  \
--sample=K2-10X5   \
--transcriptome=/mnt/mydata/refdata-gex-GRCh38-2020-A  
--localmem

#作为一款比对软件，建index肯定是必不可少的一步
#--genomeDir 索引文件夹
#--genomeFastaFiles 参考基因组fasta序列
#--sjdbGTFfile 基因注释信息gtf文件
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /Users/shijian/mydata/bulkRNA/0.reference/star.index/ --genomeFastaFiles /Users/shijian/mydata/bulkRNA/0.reference/GRCh38.p13.genome.fa --sjdbGTFfile /Users/shijian/mydata/bulkRNA/0.reference/gencode.v38.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 100

#比对
STAR --runThreadN 10 --genomeDir /Users/shijian/mydata/bulkRNA/0.reference/star.index/ --twopassMode Basic --readFilesIn /Users/shijian/mydata/bulkRNA/1.fastq/P1_1.fq.gz /Users/shijian/mydata/bulkRNA/1.fastq/P1_2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /Users/shijian/mydata/bulkRNA/3.bam/shi

#cellRanger count 
cd  /mnt/tianjin_mouse_single_cell/1.counts/mouse_J007/N
cellranger count --id=mRNA   \
--fastqs=/mnt/tianjin_mouse_single_cell/0-R848_single_cell_N_S_R_SR3/X101SC21081645-Z01-J007-sunpengya0927039/X101SC21081645-Z01-J007-B8-38_10X_release_20210924/Rawdata/N  \
--sample=N-1,N-2,N-3,N-4,N-5,N-6,N-7,N-8   \
--transcriptome=/mnt/tianjin_mouse_single_cell/0.reference_mouse/refdata-gex-mm10-2020-A  









