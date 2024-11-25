#前情回顾:
#1. hg19_to_mm10: 转换到鼠的位点(删了所有“广义保守”的行)(此步没改)
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新发现的m6A(含更多metadata)(增加了"nodsHg19_score"和"nodsmm10_score")

################################################################################
#########################      Sequence Analysis      ##########################
################################################################################
A_11bp <- DNAStringSet(Views(Hsapiens,newdata + 5)) #A_11bp("DNAStringSet"): 新发现的m6A的人11bp序列
#liftover from hg19 to mm10
chain <- import.chain('./hg19ToMm10.over.chain')
hg19_to_mm10 <- unlist(liftOver(newdata, chain)) #新的hg19_to_mm10; "CompressedGRangesList" to GRanges 

if(length(hg19_to_mm10) != 0){
  match <- match(hg19_to_mm10$ID,newdata$ID) #newdata到hg19_to_mm10有一些转不过去
  hg19_to_mm10$A_11bp <- A_11bp[match] #可以转换过去的行给上A_11bp(新发现的m6A的人11bp序列)
  hg19_to_mm10$A_mm10 <- DNAStringSet(Views(Mmusculus,hg19_to_mm10)) #鼠碱基
  
  indx <- which(hg19_to_mm10$A_mm10 == 'A') #如果甲基化位置不是“A”，则不打分
  if(length(indx) != 0){
    hg19_to_mm10 <- hg19_to_mm10[indx]
    hg19_to_mm10$A_11bp_mm10 <- DNAStringSet(Views(Mmusculus,hg19_to_mm10 + 5)) #鼠的对应位置的11bp的序列
    hg19_to_mm10$status <- '-'
    match2 <- match(hg19_to_mm10$site_pos_mm10,m6A_SB_mm10$site_pos) #匹配鼠的甲基化位点
    qhits <- which(match2 != 'NA') #hg19_to_mm10的index
    hg19_to_mm10$status[qhits] <- 'conserved' 
    
    hg19_to_mm10$score <- 0
    options(digits=2) #options()主要用于设置与计算和显示有关的全局变量,digits=2控制打印时显示的有效数字的位数为2
    
    SM_motif <- matrix(NA,4,4) #SM_motif 区域的打分矩阵
    colnames(SM_motif) <- c('A','G','C','T')
    rownames(SM_motif) <- c('A','G','C','T')
    SM_motif[1,] <- c(4,3,0,0)
    SM_motif[2,] <- c(3,4,0,0)
    SM_motif[3,] <- c(0,0,4,3)
    SM_motif[4,] <- c(0,0,3,4)
    
    SM <- matrix(NA,4,4) #其他区域的打分矩阵
    colnames(SM) <- c('A','G','C','T')
    rownames(SM) <- c('A','G','C','T')
    SM[1,] <- c(2,1,0,0)
    SM[2,] <- c(1,2,0,0)
    SM[3,] <- c(0,0,2,1)
    SM[4,] <- c(0,0,1,2)
    
    for (i in 1:length(hg19_to_mm10)) { 
      seq_hg19 <- hg19_to_mm10$A_11bp[[i]] #[[i]]取出11bp的人序列
      seq_mm10 <- hg19_to_mm10$A_11bp_mm10[[i]] #[[i]]取出11bp的鼠序列
      score_1 <- 0 
      
      for (j in 1:length(seq_hg19)) {
        #序列总长度为11bp
        if(j < 4){ #位置为-5,-4,-3
          a <- match(as.character(seq_hg19[j]),colnames(SM_motif)) 
          #match返回第一个变量中的字符匹配上第二个变量的位置, 取出human序列的碱基在SM_motif中的行号
          b <- match(as.character(seq_mm10[j]),colnames(SM_motif)) #取出mouse序列的碱基在SM_motif中的列号,colnames = rownames
          score_1 <- score_1 + SM[a,b] #加上SM对应的分数
        }
        
        if(j > 3 & j < 9){ #位置为-2,-1,0,1,2
          a <- match(as.character(seq_hg19[j]),colnames(SM_motif))
          b <- match(as.character(seq_mm10[j]),colnames(SM_motif))
          score_1 <- score_1 + SM_motif[a,b] #加上SM_motif对应的分数
        }
        
        if(j > 8){ #位置为3,4,5
          a <- match(as.character(seq_hg19[j]),colnames(SM_motif))
          b <- match(as.character(seq_mm10[j]),colnames(SM_motif))
          score_1 <- score_1 + SM[a,b]
        }
        
      }
      
      score_normalized <- score_1/32 #归一化，除以最大分数
      hg19_to_mm10$score[i] <- score_normalized
      print(i/length(hg19_to_mm10)*100) #显示打分进度
    }
    
    newdata$sequence_score <- 0
    
    match <- match(newdata$ID,hg19_to_mm10$ID)
    qhits <- which(match != 'NA')
    shits <- na.omit(match)
    newdata[qhits]$sequence_score <- hg19_to_mm10[shits]$score
    
  }else{
        newdata$sequence_score <- 0
  }
  
  
}else{
      newdata$sequence_score <- 0
}


rm(A_11bp,chain,seq_hg19,seq_mm10,SM,SM_motif)
rm(a,b,i,indx,j,match,match2,qhits,score_1,score_normalized,shits)

