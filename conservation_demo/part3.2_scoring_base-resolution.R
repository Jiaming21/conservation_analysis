#前情回顾: 
#newdata: 新发现的m6A, 含更多metadata
#site_pos_mm10: 存鼠的位点信息
#m6A_SB_mm10: 鼠保守m6A数据库

################################################################################
##########################      Base-resolution       ##########################
################################################################################
#默认0分
newdata$br_score <- 0
#保守的给1分
indx <- which(newdata$conserved_status == 'conserved')
if(length(indx) != 0){
  newdata[indx]$br_score <- 1 
}


#对鼠操作
if(length(hg19_to_mm10) != 0){
  
  #默认0分
  hg19_to_mm10$br_score <- 0 
  #删除保守行(已经打过分了) 
  match<-match(hg19_to_mm10$site_pos_mm10,m6A_SB_mm10$site_pos) #鼠信息之间相互比较
  qhits <- which(match != 'NA')  
  if(length(qhits) != 0){
    hg19_to_mm10 <- hg19_to_mm10[-qhits] 
  }
  
  #给偏移符合的打分
  list <- c(0.8,0.6,0.4,0.2)
  for (i in 1:4) {
    overlap <- findOverlaps(hg19_to_mm10 + i,m6A_SB_mm10) 
    #overlap: 输出偏移x基因位置匹配保守库的鼠行与保守行；queryHits：对应鼠的行号，subjectHits：对应鼠保守数据库的行号
    if(length(overlap) != 0){
      hg19_to_mm10[queryHits(overlap)]$br_score <- list[i]  
      #把信息同步到newdata上
      match <- match(hg19_to_mm10[queryHits(overlap)]$ID, newdata$ID) 
      shits <- na.omit(match)
      newdata$br_score[shits] <- list[i]
      newdata$PubmedID_mm10[shits] <- m6A_SB_mm10[subjectHits(overlap)]$PubmedID
      newdata$GSE_mm10[shits] <- m6A_SB_mm10[subjectHits(overlap)]$GSE
      newdata$Technique_mm10[shits] <- m6A_SB_mm10[subjectHits(overlap)]$Technique
      newdata$Cell_line_mm10[shits] <- m6A_SB_mm10[subjectHits(overlap)]$Cell_line
      newdata$Treatment_mm10[shits] <- m6A_SB_mm10[subjectHits(overlap)]$Treatment
      newdata$NumberOfExperimentSupported_mm10[shits] <- m6A_SB_mm10[subjectHits(overlap)]$NumberOfExperimentSupported
      newdata$ID_mm10[shits] <- m6A_SB_mm10[subjectHits(overlap)]$ID
      #从hg19_to_mm10中移除匹配上的位点
      hg19_to_mm10 <- hg19_to_mm10[-queryHits(overlap)]
    }
    print(i)
    
  }
  
  
}


rm(overlap)
rm(i)
rm(indx)
rm(list)
rm(match)
rm(qhits)
rm(shits)
