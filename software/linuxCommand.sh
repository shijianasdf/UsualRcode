#
# 2014.8
#
#
# 此文件用于系统初始化，运行完毕后，软件方面将会恢复
# 另外此时仍需对firefox安装adobe flash player插件
# 需要事先把一些软件复制到Downloads/Software中
# 使用root权限执行
#
#
#
# 需要固定的目录储存一些软件包和配置文件
# 需要存在 $sdaDestDir/Software 目录下
# /sda/Software/mysql.my.cnf 也要存在


##
## 配置软件源  ----  开始
##

#yum-config-manager --add-repo=https://copr.fedoraproject.org/coprs/mosquito/myrepo/repo/epel-$(rpm -E %?rhel)/mosquito-myrepo-epel-$(rpm -E %?rhel).repo
yum install -y epel-release.noarch
# 因为centos7没有对应的rpmfushion，centos7是不能安裝el6的源的
# 因此这里会有小问题，导致后面的copr安装不上
# 但copr官网依旧指定如此安装，奇怪，未解
#yum localinstall http://li.nux.ro/download/nux/dextop/el$(rpm -E %rhel)/x86_64/nux-dextop-release-0-2.el$(rpm -E %rhel).nux.noarch.rpm http://download1.rpmfusion.org/nonfree/el/updates/$(rpm -E %rhel)/x86_64/rpmfusion-nonfree-release-$(rpm -E %rhel)-1.noarch.rpm http://download1.rpmfusion.org/free/el/updates/$(rpm -E %rhel)/x86_64/rpmfusion-free-release-$(rpm -E %rhel)-1.noarch.rpm 

# 注意：rhel/centos 6/7 需要添加 epel 、 rpmfusion 、Nux Dextop 源。Nux 源的个别包与 base 源有冲突，建议使用 yum-plugin-priorities 为源分级。 
# 但是没有弄。。。
#yum install -y yum-plugin-priorities
#vim /etc/yum.repos.d/RepoName.repo
#[repo_name]
#name= Repo Full Name
#baseurl= Repo URL
#enabled= 1(enable) 0(disable)
#priority= priority number (range: 1-99, 1 high priority)
# 建议 base, updates, epel, mosquito-myrepo 优先级为 1，其他源 (rpmfusion, remi, Nux Dextop, RPMforge 等) 设为 2。这样会减少源之间的软件包冲突。
# 更新时有冲突的，可使用 # yum update --exclude=Package_Name 来排除某个软件包。
# 因安装脚本报错而无法卸载的，可使用 # rpm -e --noscripts Package_Name 来卸载软件包，目前正在除虫。

# yum install -y dnf
# 这个是一个新的软件包管理软件，类似于yum，主要是为了调用copr
# 其目的在于安装翻墙软件shadowsocks，但没有安装上。。。
# 因此注释掉
# 另
# https://copr.fedoraproject.org/coprs/mosquito/myrepo/
# 这个copr软件源中有好东西

##
## 配置软件源  ----  完毕
##



##
## yum安装软件  ----  开始
##

# NTFS支持
yum install -y ntfs-3g ntfsprogs
# 为了编译安装mysql
yum install -y bison-devel bison cmake gcc-c++ gcc gcc-gfortran ncurses-devel unrar
# 为了安装虚拟机VirtualBox
yum install -y kernel-devel kernel
# 基础软件
yum install -y httpd gimp tigervnc nmap kate perl git R R-XML libcurl-devel vim gtk2 gtk2-devel iftop nmap iptraf k3b filezilla vsftpd vsftpd-sysvinit
# 下载工具 (不好用。。。)
# yum install -y pointdownload
# 百度盘 (不好用。。。)
# yum install -y bcloud
# 音频和视频播放器
yum install -y vlc vlc-extras audacious audacious-plugins
# 为了翻墙
# 因为中共的网络封锁可能会出现网络连线的错误
wget -i https://copr.fedoraproject.org/coprs/librehat/shadowsocks/repo/epel-7/librehat-shadowsocks-epel-7.repo  在哪个目录就下载到哪个目录
cp librehat-shadowsocks-epel-7.repo /etc/yum.repos.d/
rm -f librehat-shadowsocks-epel-7.repo
yum install -y shadowsocks-qt5

##
## yum安装软件  ----  完毕
##



##
## rpm安装软件  ----  开始
##
# 重要！！！
# 需要把Software事先复制到Home目录的Downloads中
cd Downloads/Software/
yum install -y VirtualBox-5.0-5.0.0_101573_el7-1.x86_64.rpm
yum install -y google-chrome-stable_current_x86_64.rpm
##
## rpm安装软件  ----  开始
##



##
## source安装软件  ----  开始
##


# 设定临时的环境变量，并建立目录    ---- 开始
# 重要！！！
# 如若想改变，改此即可
export customDestDir=/home/sj
# 重要！！！
# 备份profile文件和PATH
cp /etc/profile /etc/profile.bk
echo "export PATH=$PATH" >> /etc/profile.bk
echo 'export customDestDir=/home/sj' >> /etc/profile

# 将建立目录
mkdir $customDestDir/software
mkdir $customDestDir/mysqlData
# 设定临时的环境变量，并建立目录    ---- 完毕

## Software 目录安装    ----  开始
# 安装 TOMCAT 7
tar zxvf apache-tomcat-7.0.55.tar.gz
mv apache-tomcat-7.0.55 $customDestDir/Software/
echo 'export CATALINAHOME=$customDestDir/Software/apache-tomcat-7.0.55' >> /etc/profile

# 安装 ECLIPSE
tar zxvf eclipse-jee-mars-1-linux-gtk-x86_64.tar.gz
mv eclipse $customDestDir/Software/

# 安装 JDK 1.7
tar zxvf jdk-7u65-linux-x64.tar.gz
mv jdk1.7.0_65 $customDestDir/Software/
echo 'export JAVAHOME=$customDestDir/Software/jdk1.7.0_65' >> /etc/profile

# 安装 ASPERA
#tar zxvf aspera-connect-3.5.1.92523-linux-64.tar.gz 
#chmod +x aspera-connect-3.5.1.92523-linux-64.sh 
#./aspera-connect-3.5.1.92523-linux-64.sh 
#mv /root/.aspera $customDestDir/Software/
#mv $customDestDir/Software/.aspera $customDestDir/Software/aspera
#echo 'export ASPERAHOME=$customDestDir/Software/aspera/connect' >> /etc/profile

## Software 目录安装    ----  结束

## BioSoftware 目录安装    ----  开始
#cd $customDestDir/BioSoftware
# BEDTOOLS
#git clone https://github.com/arq5x/bedtools2.git
#cd bedtools2/
#make
#echo 'export BEDTOOLS2HOME=$customDestDir/BioSoftware/bedtools2' >> /etc/profile

# SAMTOOLS
#git clone https://github.com/samtools/samtools.git
# 新版本的samtools需要这个
# 并没有编译这个？？？
#git clone https://github.com/samtools/htslib.git
#cd samtools
#make
#echo 'export SAMTOOLS=$customDestDir/BioSoftware/samtools' >> /etc/profile
# circos-0.67
# tophat2
# UCSCBin
# bowtie
# Cufflinks
# bowtie2       
# FastQC
# ncbi-blast-2.2.30+
# sratoolkit




## BioSoftware 目录安装    ----  结束



# 向profile中添加各软件包的环境变量
# 这里面的MYSQLHOME是为了mysql的，在随后的编译中将产生相关的文件目录
# 注意此时MYSQLHOME还没有，只是把它添加了而已
echo 'export MYSQLHOME=/usr/local/mysql' >> /etc/profile
#echo 'export PATH=$SAMTOOLS:$BEDTOOLS2HOME/bin:$ASPERAHOME/bin:$JAVAHOME/bin:$CATALINAHOME/bin:$MYSQLHOME/bin:$PATH' >> /etc/profile
echo 'export PATH=$JAVAHOME/bin:$CATALINAHOME/bin:$MYSQLHOME/bin:$PATH' >> /etc/profile
echo 'export CLASSPATH=.:$JAVAHOME/lib/' >> /etc/profile

# 修改所有者和用户组
chown -R sda $customDestDir
chgrp -R sda $customDestDir

# mysql编译安装 需要交互，最后安装这个
# 上面安装过，所以注释了
# yum -y install bison-devel bison cmake gcc-c++ gcc gcc-gfortran ncurses-devel
tar zxvf mysql-5.6.20.tar.gz
cd mysql-5.6.20
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/mysql -DMYSQL_DATADIR=/sj/mysqlData/ -DSYSCONFDIR=/etc -DENABLED_LOCAL_INFILE=1 -DWITH_INNOBASE_STORAGE_ENGINE=1 -DWITH_ARCHIVE_STORAGE_ENGINE=1 -DWITH_BLACKHOLE_STORAGE_ENGINE=1 -DWITH_EXAMPLE_STORAGE_ENGINE=1 -DWITH_FEDERATED_STORAGE_ENGINE=1 -DWITH_INNOBASE_STORAGE_ENGINE=1 -DWITH_PARTITION_STORAGE_ENGINE=1 -DWITH_PERFSCHEMA_STORAGE_ENGINE=1 -DWITH_READLINE=1 -DEXTRA_CHARSETS=all -DDEFAULT_CHARSET=utf8 -DDEFAULT_COLLATION=utf8_general_ci -DMYSQL_UNIX_ADDR=/sj/mysqlData/mysql.sock
make
make install
cd /usr/local/mysql
groupadd mysql
useradd -r -g mysql mysql
chown -R mysql .
chgrp -R mysql .
scripts/mysql_install_db --basedir=/usr/local/mysql --datadir=/sj/mysqlData --user=mysql --ldata=/sj/mysqlData
chown -R root .
chown -R mysql /sj/mysqlData
mkdir /var/log/mariadb
bin/mysqld_safe --user=mysql &
cp support-files/mysql.server /etc/init.d/mysqld
cp $customDestDir/Software/mysql.my.cnf /etc/my.cnf
bin/mysqladmin -uroot -p password
##
## source安装软件  ----  完毕
##


#翻墙
yum install -y openvpn
cd /sj/Software/openvpn
openvpn --config localhost.localdomain.ovpn  


#修改目录名字为英文
export LANG=en_US
xdg-user-dirs-gtk-update

#mplayer安装
CentOS官方源里的电影播放器是根本就用不了的。所以这里我们要利用第三方源来安装。加入RPMforge源：
sudo rpm -Uhv http://rpmforge.sw.be/redhat/el5/en/i386/rpmforge/RPMS/rpmforge-release-0.3.6-1.el5.rf.i386.rpm 
yum install mplayer kmplayer
编辑/etc/yum.repos.d/rpmforge.repo文件，将其中的enabled项的值设成1表示启用源：

#启动mysql
根用户：service mysql start
普通用户：mysql -uroot -p123456

#Linux下一个非常强大的多线程下载工具aria2c
yum search aria
yum install aria2.x86_64
aria2c downloadFiles 
aria2c -i downloadFiles

#aria2c用法
aria2c http://AAA.BBB.CCC/file.zip   普通下载
aria2c -s 2 http://AAA.BBB.CCC/file.zip  开2个线程下载
aria2c http://AAA.BBB.CCC/file.zip ftp://DDD.EEE.FFF/GGG/file.zip  从不同的地址下载同一文件
aria2c http://AAA.BBB.CCC/file.zip ftp://DDD.EEE.FFF/GGG/file.zip  支持不同的协议下载同一文件
aria2c -o test.torrent http://AAA.BBB.CCC/file.torrent  下载BT种子
aria2c --max-upload-limit 40K -T file.torrent  设定BT最大上传速度
aria2c http://AAA.BBB.CCC/file.metalink   从metalink下载文件

#软件源安装好处
相比于源码编译安装，软件源yum安装直接自动下载安装好dependencies，并且yum安装自动设置好环境变量

#linux文件名乱码解决方法
文件是在Windows下创建的，而Windows的文件名中文编码默认GBK，Linux中默认文件名编码为UTF-8，编码不一致导致了文件名乱码的问题，解决这个问题需要对文件名进行转码，这个工具就是convmv。
convmv软件下载: yum install -y convmv.noarch
convmv软件使用方法：convmv -f GBK -t zh_CN.UTF-8 -r --notest ./*   将当前目录下的所有文件从GBK格式改为UTF-8


#查看磁盘空间
df -h

