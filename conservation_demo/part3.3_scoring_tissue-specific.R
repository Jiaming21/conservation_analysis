#前情回顾: 
#newdata: 新发现的m6A(含更多metadata)

################################################################################
##########################      Tissue specific      ###########################
################################################################################

newdata$brts_score <- 0

list <- c("brain","kidney","liver","ESC")

for (i in 1:4) {
  list_hg19 <- grep(list[i],newdata$Cell_line) #grep：查找目标字符串或字符串向量中是否包含目标子串,返回位置信息index
  list_mm10 <- grep(list[i],newdata$Cell_line_mm10)
  a <- intersect(list_hg19,list_mm10) #intersect: 两个向量的交集，集合可以是数字、字符串等,返回对应的值
  
  if(length(a) != 0){
    newdata[a]$brts_score <- newdata[a]$br_score #对应上给分与br_score一样
  }
  
}

rm(a,i,list,list_hg19,list_mm10)


