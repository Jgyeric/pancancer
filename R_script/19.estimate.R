#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#引用包
library(limma)
library(estimate)
setwd("D:\\Rdata\\panCancer\\19.estimate")        #设置工作目录

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
  out=data[rowMeans(data)>0,]
  out=rbind(ID=colnames(out),out)
  #输出整理后的矩阵文件
  write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)
  
  #运行estimate包
  filterCommonGenes(input.f="uniq.symbol.txt", output.f="commonGenes.gct", id="GeneSymbol")
  estimateScore(input.ds = "commonGenes.gct",output.ds="estimateScore.gct")
  
  #输出每个样品的打分
  scores=read.table("estimateScore.gct",skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  rownames(scores)=gsub("\\.","\\-",rownames(scores))
  outTab=rbind(outTab,cbind(scores,CancerType))
  file.remove("commonGenes.gct")
  file.remove("estimateScore.gct")
  file.remove("uniq.symbol.txt")
}
out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="estimateScores.txt",sep="\t",quote=F,row.names=F)