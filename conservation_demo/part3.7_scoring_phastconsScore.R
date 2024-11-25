################################################################################
##########################      PhastCons score      ###########################
################################################################################

library(phastCons100way.UCSC.hg19) #The phastCons 100-way conservation scores were calculated for human genome derived from genome-wide multiple alignments with 99 other vertebrate species.
gsco <- phastCons100way.UCSC.hg19

newdata <- gscores(gsco, newdata,pop="DP2") #计算指定区域的平均保守性得分, pop 选择一个population

indx_na <- as.numeric(attr(na.omit(newdata$DP2),"na.action")) #attr(na.omit(),"na.action"):返回向量a中元素为NA的下标

if(length(indx_na) != 0){
  
  newdata$DP2[indx_na] <- 0 #将NA赋值为0
  
}else{
  
  newdata <- readRDS('/Users/huangjiaming/Desktop/conservation/webserver/example_newdata.rds')
  
}











