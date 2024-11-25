#前情回顾: 
#newdata: 新发现的m6A, 只有位置信息, 无metadata

################################################################################
###########################      Liftover      #################################
################################################################################

if(length(newdata) != 0){ #存在新发现的m6A
  newdata$PubmedID <- '-' 
  newdata$GSE <- '-'
  newdata$Technique <- '-'
  newdata$Cell_line <- '-'
  newdata$Treatment <- '-'
  newdata$site_pos <- paste0(seqnames(newdata),':',start(newdata),' ',strand(newdata)) #拼接位置信息 
  newdata$reference_sequence <- DNAStringSet(Views(Hsapiens,newdata+20)) 
  newdata$A_41bp <- ranges(Views(Hsapiens,newdata + 20))
  newdata$NumberOfExperimentSupported <- 1
  newdata$ID <- paste0('Human_newly_identified_m6A_',1:length(newdata))
  
  chain <- import.chain('./hg19ToMm10.over.chain')
  hg19_to_mm10 <- unlist(liftOver(newdata, chain)) #经转换得到的鼠的位点,metadata都是人的(人鼠一样)
  
  if(length(hg19_to_mm10) != 0){ #人中的位点有可以转换到鼠的
    
    #找到人中哪些行匹配上了鼠中的哪些行, 对匹配上的人中的行: liftoverToMm10是YES', site_pos_mm10填好，其他行'-'
    match <- match(newdata$site_pos,hg19_to_mm10$site_pos) #match=找匹配行
    qhits <- which(match != 'NA')  
    shits <- na.omit(match) 
    newdata$liftoverToMm10 <- 'NO' 
    newdata$liftoverToMm10[qhits] <- 'YES'
    newdata$site_pos_mm10 <- '-' 
    hg19_to_mm10$site_pos_mm10 <- paste0(seqnames(hg19_to_mm10),":",start(hg19_to_mm10)," ",strand(hg19_to_mm10)) 
    newdata$site_pos_mm10[qhits] <- hg19_to_mm10$site_pos_mm10[shits] 
    
    
    
    #找哪些转换成功的行又在鼠保守库里有信息
    m6A_SB_mm10 <- readRDS('./m6A_SB_mm10.rds') #鼠保守m6A数据库
    match <- match(newdata$site_pos_mm10,m6A_SB_mm10$site_pos) 
    qhits <- which(match != 'NA')
    shits <- na.omit(match)
    newdata$conserved_status <- 'not conserved' #默认全部都不保守，只有可以转换，又在鼠保守库里面才保守
    if(length(qhits) != 0){ #有可以转换，又在鼠保守库里面的
      newdata$conserved_status[qhits] <- 'conserved' #符合可以转换，又在鼠保守库里面
      
      newdata$PubmedID_mm10 <- '-' 
      newdata$PubmedID_mm10[qhits] <- m6A_SB_mm10$PubmedID[shits] 
      
      newdata$GSE_mm10 <- '-'
      newdata$GSE_mm10[qhits] <- m6A_SB_mm10$GSE[shits]
      
      newdata$Technique_mm10 <- '-'
      newdata$Technique_mm10[qhits] <- m6A_SB_mm10$Technique[shits]
      
      newdata$Cell_line_mm10 <- '-'
      newdata$Cell_line_mm10[qhits] <- m6A_SB_mm10$Cell_line[shits]
      
      newdata$Treatment_mm10 <- '-'
      newdata$Treatment_mm10[qhits] <- m6A_SB_mm10$Treatment[shits]
      
      newdata$NumberOfExperimentSupported_mm10 <- '-'
      newdata$NumberOfExperimentSupported_mm10[qhits] <- m6A_SB_mm10$NumberOfExperimentSupported[shits]
      
      newdata$ID_mm10 <- '-'
      newdata$ID_mm10[qhits] <- m6A_SB_mm10$ID[shits]
    }else{ #可以转换，但不在鼠保守库
      newdata$PubmedID_mm10 <- '-'
      newdata$GSE_mm10 <- '-'
      newdata$Technique_mm10 <- '-'
      newdata$Cell_line_mm10 <- '-'
      newdata$Treatment_mm10 <- '-'
      newdata$NumberOfExperimentSupported_mm10 <- '-'
      newdata$ID_mm10 <- '-'
    }
    
    
  }else{ #都不能转换到鼠
    newdata$liftoverToMm10 <- 'NO'
    newdata$site_pos_mm10 <- '-'
    newdata$conserved_status <- 'not conserved'
    newdata$PubmedID_mm10 <- '-'
    newdata$GSE_mm10 <- '-'
    newdata$Technique_mm10 <- '-'
    newdata$Cell_line_mm10 <- '-'
    newdata$Treatment_mm10 <- '-'
    newdata$NumberOfExperimentSupported_mm10 <- '-'
    newdata$ID_mm10 <- '-'
  }
 
   
}



rm(chain)
rm(match)
rm(qhits)
rm(shits)



