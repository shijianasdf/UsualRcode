#P3_651038下载
{
cd Users/biofly/project/shijian/macnd
./macnd
./macnd login -u X101SC21081645-Z01-J010 -p password
./macnd list
./macnd list oss://CP2021081200153
#批量下载
cat /Users/biofly/project/shijian/macnd/X101SC21081645-Z01-J010.txt | while read it
do
  command=./macnd cp oss://CP2021081200153/${it} /Volumes/HDD3/tianjin_ICC/P3_651038/fastq/X101SC21081645-Z01-J010
  echo $command >> X101SC21081645-Z01-J010.sh
done
#批量检查
ls /Volumes/HDD3/tianjin_ICC/P3_651038/fastq/X101SC21081645-Z01-J010/*.fastq.gz | while read id
do
  md5 $id
done
#批量下载
./macnd login -u X101SC21081645-Z01-J011 -p password
./macnd list
./macnd list oss://CP2021081200153 > X101SC21081645-Z01-J011.txt(手动copy过去)
cat /Users/biofly/project/shijian/macnd/X101SC21081645-Z01-J011.txt | while read it
do
  command=./macnd cp oss://CP2021081200153/${it} /Volumes/HDD3/tianjin_ICC/P3_651038/fastq/X101SC21081645-Z01-J011
  echo $command >> X101SC21081645-Z01-J011.sh
done
#批量检查
ls /Volumes/HDD3/tianjin_ICC/P3_651038/fastq/X101SC21081645-Z01-J011/*.fastq.gz | while read id
do
  md5 $id
done
#批量下载
./macnd login -u X101SC21081645-Z01-J012 -p password
./macnd list
./macnd list oss://CP2021081200153 > X101SC21081645-Z01-J012.txt(手动copy过去)
cat /Users/biofly/project/shijian/macnd/X101SC21081645-Z01-J012.txt | while read it
do
  command=./macnd cp oss://CP2021081200153/${it} /Volumes/HDD3/tianjin_ICC/P3_651038/fastq/X101SC21081645-Z01-J012
  echo $command >> X101SC21081645-Z01-J012.sh
done
#批量检查
ls /Volumes/HDD3/tianjin_ICC/P3_651038/fastq/X101SC21081645-Z01-J012/*.fastq.gz | while read id
do
  md5 $id
done
