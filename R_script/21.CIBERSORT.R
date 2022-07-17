######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#引用包
library(limma)
setwd("D:\\Rdata\\panCancer\\21.CIBERSORT")        #设置工作目录
pFilter=0.05

#读取目录下的文件
files=dir()
files=grep("^symbol.",files,value=T)

outTab=data.frame()
for(i in files){
	#读取文件
	CancerType=gsub("symbol\\.|\\.txt","",i)
	rt=read.table(i,sep="\t",header=T,check.names=F)
	
	#如果一个基因占了多行，取均值
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	
	#删除正常，只保留肿瘤样品
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2","1",group)
	data=data[,group==0]
	data=data[rowMeans(data)>0,]

	#数据矫正
	v <-voom(data, plot = F, save.plot = F)
	out=v$E
	out=rbind(ID=colnames(out),out)
	write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #输出文件
	
	#运行CIBERSORT，得到免疫细胞含量结果
	source("panCancer21.CIBERSORT.R")
	results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)

	#输出每个样品的打分
	immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
	immune=immune[immune[,"P-value"]<pFilter,]
	immune=as.matrix(immune[,1:(ncol(immune)-3)])
	outTab=rbind(outTab,cbind(immune,CancerType))
	file.remove("CIBERSORT-Results.txt")
	file.remove("uniq.symbol.txt")
}
out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="CIBERSORT.result.txt",sep="\t",quote=F,row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
