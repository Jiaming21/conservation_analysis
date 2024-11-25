newdata$tissue_specific <- '-'

tissueList <- c("brain","kidney","liver","ESC" )

for (i in 1:4) {
  
  brain_list_hg19 <- grep(tissueList[i],newdata$Cell_line) #grep：查找目标字符串或字符串向量中是否包含目标子串,返回位置信息index
  brain_list_mm10 <- grep(tissueList[i],newdata$Cell_line_mm10)
  a <- intersect(brain_list_hg19,brain_list_mm10) #intersect:两个向量的交集，集合可以是数字、字符串等,返回对应的值
  newdata$tissue_specific[a] <- 'YES'
  
}