##Rstudio远程连接服务器linux
ssh shijian2015@210.46.80.170 (Im1)
ssh shijian2015@210.46.80.168:1 (XteamServer69!~)
/pub5/xiaoyun/Bin/R
##xshell远程登录服务器linux
xteam
ip
c1(111111)
R2
#进入论坛函数
library(RSc)
RS("b1")
test <- RSE(IP("IPD.queryGenomics.ser","UCSC.CODING.GENE.gene",genome="hg19"));
head(test)
# 退出
RSC()
# 重启b1节点
RSS("b1")
RSS("x1")
##远程连接服务器源码
git@210.46.80.150:/pub5/xiaoyun/Software/IP.GIT
xteam104


# 备份自己的环境变量
mv ~/.bash_profile ~/.bash_profile.bk
mv ~/.bashrc ~/.bashrc.bk
## ~ 家目录 将刘威的bash_profile和bashrc复制到自己的家目录中
cp /pub6/Temp/liuwei/.bash_profile ~/.bash_profile
cp /pub6/Temp/liuwei/.bashrc ~/.bashrc
