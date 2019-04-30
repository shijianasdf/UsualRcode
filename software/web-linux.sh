#tigerVNC服务器地址和密码
202.97.205.69:2
XiaHMU~217LXZ1644!han~IP

2018.05.28-new
host:202.97.205.69:1
name：admin_s
password：XteamServer69!~

#删除已有的数据库
sudo rm -r daec

#修改daec文件夹及所有子文件的所有者和用户组全为mysql
sudo chown -R mysql:mysql /pub5/ftp/public/LnChrom/DEC/daec

#copy daec到mysqldata
sudo cp -r /pub5/ftp/public/LnChrom/DEC/daec  /pub5/xiaoyun/MySQLData 

#creat usr
mysql -uroot -ppassWordis204!
GRANT all privileges on daec.* To daec@localhost IDENTIFIED BY 'daec_xteam'; #用户名:daec 密码:daec_xteam
exit

#修改MySQL的参数 my.cnf文件来解决Mysql连接超时的问题，wait_timeout最大为31536000即1年，在my.cnf中加入：
[mysqld]
wait_timeout=31536000
interactive_timeout=31536000
重启生效，需要同时修改这两个参数。
show variables like '%timeout%';

#重启Mysql服务
sudo service mysql.server restart 
sudo service mysqld stop
sudo service mysqld start 




