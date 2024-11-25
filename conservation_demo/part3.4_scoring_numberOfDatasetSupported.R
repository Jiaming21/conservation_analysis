#前情回顾: newdata: 新发现的m6A, 含更多的metadata

################################################################################
####################      Number of dataset supported      #####################
################################################################################
#hg19
newdata$nodsHg19_score <- 0
#基于当前研究得到的结果--赋值为0.2
newdata[which(newdata$NumberOfExperimentSupported == 1)]$nodsHg19_score <- 0.2

#mm10
newdata$nodsmm10_score <- 0
#which的用法, 返回index 信息
if(length(which(newdata$NumberOfExperimentSupported_mm10 == 1)) != 0){
  newdata[which(newdata$NumberOfExperimentSupported_mm10 == 1)]$nodsmm10_score <- 0.2
}

if(length(which(newdata$NumberOfExperimentSupported_mm10 == 2)) != 0){
  newdata[which(newdata$NumberOfExperimentSupported_mm10 == 2)]$nodsmm10_score <- 0.4
}

if(length(which(newdata$NumberOfExperimentSupported_mm10 >= 3 & newdata$NumberOfExperimentSupported_mm10 <= 4)) != 0){
  newdata[which(newdata$NumberOfExperimentSupported_mm10 >= 3 & newdata$NumberOfExperimentSupported_mm10 <= 4)]$nodsmm10_score <- 0.6
}

if(length(which(newdata$NumberOfExperimentSupported_mm10 >= 5 & newdata$NumberOfExperimentSupported_mm10 <= 6)) != 0){
  newdata[which(newdata$NumberOfExperimentSupported_mm10 >= 5 & newdata$NumberOfExperimentSupported_mm10 <= 6)]$nodsmm10_score <- 0.8
}

if(length(which(newdata$NumberOfExperimentSupported_mm10 > 6)) != 0){
  newdata[which(newdata$NumberOfExperimentSupported_mm10 > 6)]$nodsmm10_score <- 1
}
