################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/score/hg38ToMm10/file_A_hg38ToMm10_sequenceAnalysis')
################################################################################
#########################      Sequence Analysis      ##########################
################################################################################
hg38_tombo_0109 <- readRDS('hg38_tombo_0109.rds')
hg38_tombo_0109 <- hg38_tombo_0109[which(vcountPattern('single base',hg38_tombo_0109$resolution)>=1)]
hg38_tombo_0109 <- hg38_tombo_0109[which(hg38_tombo_0109$ref_base=='A')]
hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))

chain <- import.chain('hg38ToMm10.over.chain')
hg38_to_mm10 <- unlist(liftOver(hg38_tombo_0109,chain))
hg38_to_mm10$site_pos_mm10 <- paste0(seqnames(hg38_to_mm10),':',start(hg38_to_mm10),' ',strand(hg38_to_mm10))

match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$liftOverToMm10 <- 'No'
hg38_tombo_0109$liftOverToMm10[qhits] <- 'Yes'

mm10_tombo_0109 <- readRDS('./mm10_tombo_0109.rds')
mm10_tombo_0109 <- mm10_tombo_0109[which(vcountPattern('single base',mm10_tombo_0109$resolution)>=1)]
mm10_tombo_0109 <- mm10_tombo_0109[which(mm10_tombo_0109$ref_base=='A')]
mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))

match <- match(hg38_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
qhits <- which(match != 'NA')
shits <- na.omit(match)
hg38_to_mm10$status <- 'not conserved'
hg38_to_mm10$status[qhits] <- 'conserved'

hg38_to_mm10 <- hg38_to_mm10[hg38_to_mm10$status == 'conserved']

hg38_tombo_0109$A_31bp_hg38 <- DNAStringSet(Views(Hsapiens,hg38_tombo_0109 + 15))
match <- match(hg38_to_mm10$site_pos_hg38,hg38_tombo_0109$site_pos_hg38)
qhits <- which(match != 'NA')
shits <- na.omit(match)
hg38_to_mm10$A_31bp_hg38 <- DNAStringSet('-')
hg38_to_mm10$A_31bp_hg38[qhits] <- hg38_tombo_0109$A_31bp_hg38[shits]
hg38_to_mm10$A_31bp_mm10 <- DNAStringSet(Views(Mmusculus,hg38_to_mm10 + 15))

#-------------------------------------------------------------------------------
hg38_to_mm10$score <- 0

SM_1 <- matrix(NA,4,4) #SM_1 区域的打分矩阵
colnames(SM_1) <- c('A','G','C','T')
rownames(SM_1) <- c('A','G','C','T')
SM_1[1,] <- c(5,4,0,0)
SM_1[2,] <- c(4,5,0,0)
SM_1[3,] <- c(0,0,5,4)
SM_1[4,] <- c(0,0,4,5)

SM_2 <- matrix(NA,4,4) #SM_2 区域的打分矩阵
colnames(SM_2) <- c('A','G','C','T')
rownames(SM_2) <- c('A','G','C','T')
SM_2[1,] <- c(4,3,0,0)
SM_2[2,] <- c(3,4,0,0)
SM_2[3,] <- c(0,0,4,3)
SM_2[4,] <- c(0,0,3,4)

SM_3 <- matrix(NA,4,4) #打分矩阵
colnames(SM_3) <- c('A','G','C','T')
rownames(SM_3) <- c('A','G','C','T')
SM_3[1,] <- c(2,1,0,0)
SM_3[2,] <- c(1,2,0,0)
SM_3[3,] <- c(0,0,2,1)
SM_3[4,] <- c(0,0,1,2)

for (i in 1:length(hg38_to_mm10)) {
  seq_hg38 <- hg38_to_mm10$A_31bp_hg38[[i]]
  seq_mm10 <- hg38_to_mm10$A_31bp_mm10[[i]]
  score_1 <- 0 
  
  for (j in 1:length(seq_hg38)) {
    
    if(j < 6){ 
      a <- match(as.character(seq_hg38[j]),rownames(SM_3)) 
      b <- match(as.character(seq_mm10[j]),colnames(SM_3)) 
      score_1 <- score_1 + SM_3[a,b] 
    }
    
    if(j > 5 & j < 11){
      a <- match(as.character(seq_hg38[j]),rownames(SM_2)) 
      b <- match(as.character(seq_mm10[j]),colnames(SM_2)) 
      score_1 <- score_1 + SM_2[a,b] 
    }
    
    if(j > 10 & j < 22){
      a <- match(as.character(seq_hg38[j]),colnames(SM_1))
      b <- match(as.character(seq_mm10[j]),colnames(SM_1))
      score_1 <- score_1 + SM_1[a,b]
    }
    
    if(j > 21 & j < 27){
      a <- match(as.character(seq_hg38[j]),colnames(SM_2))
      b <- match(as.character(seq_mm10[j]),colnames(SM_2))
      score_1 <- score_1 + SM_2[a,b]
    }
    
    if(j > 26){
      a <- match(as.character(seq_hg38[j]),colnames(SM_3))
      b <- match(as.character(seq_mm10[j]),colnames(SM_3))
      score_1 <- score_1 + SM_3[a,b]
    }
    
  }
  
  score_normalized <- score_1/125 #归一化，除以最大分数
  hg38_to_mm10$sequence_score[i] <- score_normalized
  print(i/length(hg38_to_mm10)*100) #显示打分进度
  
}

match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA')
shits <- na.omit(match)
hg38_tombo_0109$sequence_score <- 0
hg38_tombo_0109[qhits]$sequence_score <- hg38_to_mm10[shits]$sequence_score

#-------------------------------------------------------------------------------
saveRDS(hg38_tombo_0109,'../file/hg38_tombo_0109_sa.rds')



