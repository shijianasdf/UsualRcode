#目前存在大量的注释信息的数据库，我们需要一个方便的搜索工具，用于找到我们所需要的信息。Biconductor建立在R语言上的一个开源项目，旨在未高通量数据分析提供可靠的工具。项目的一个重要部分就是组织网络上大量的注释信息，方便科研人员使用。
#目前最新的工具包叫做AnnotationHub，顾名思义，就是注释信息的中装站。通过它，能找到了几乎所有的注释资源。如果没有，你还可以根据已有的数据用它提供的函数进行构建。
library(AnnotationHub)
ah <- AnnotationHub()
ah$dataprovider #查看数据来源
ah$species
ah$rdataclass
mcols(ah) #更多的元信息
ath <- ah[['AH53758']] #获取数据
query(ah, c("OrgDb","Homo sapiens"), ignore.case=TRUE)