################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/score/hg38ToMm10/file_A_hg38ToMm10_base-resolution')
################################################################################
############################      base-resolution    ###########################
################################################################################
hg38_tombo_0109 <- readRDS('hg38_tombo_0109.rds')
hg38_tombo_0109 <- hg38_tombo_0109[which(vcountPattern('single base',hg38_tombo_0109$resolution)>=1)]
hg38_tombo_0109 <- hg38_tombo_0109[which(hg38_tombo_0109$ref_base=='A')]
hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))

chain <- import.chain('hg38ToMm10.over.chain')
hg38_to_mm10 <- unlist(liftOver(hg38_tombo_0109,chain))
match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$liftOverToMm10 <- 'No'
hg38_tombo_0109$liftOverToMm10[qhits] <- 'Yes'

mm10_tombo_0109 <- readRDS('mm10_tombo_0109.rds')
mm10_tombo_0109 <- mm10_tombo_0109[which(vcountPattern('single base',mm10_tombo_0109$resolution)>=1)]
hg38_to_mm10$site_pos_mm10 <- paste0(seqnames(hg38_to_mm10),':',start(hg38_to_mm10),' ',strand(hg38_to_mm10))
mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))

match <- match(hg38_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_to_mm10$modification_mm10 <- '-'
hg38_to_mm10$modification_mm10[qhits] <- mm10_tombo_0109$ref_base[shits]

#-------------------------------------------------------------------------------

match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10[hg38_to_mm10$modification_mm10 == 'A']$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$br_score <- 0
if(length(qhits) != 0){
  hg38_tombo_0109[qhits]$br_score <- 1 
}


if(length(hg38_to_mm10) != 0){
  
  hg38_to_mm10$br_score <- 0 

  if(length(qhits) != 0){
    hg38_to_mm10 <- hg38_to_mm10[-shits] #删除保守行(已经打过分了) 
  }
  
  #给偏移符合的打分
  list <- c(0.8,0.6,0.4,0.2)
  for (i in 1:4) {
    overlap <- findOverlaps(hg38_to_mm10 + i,mm10_tombo_0109[mm10_tombo_0109$ref_base == 'A'])
    if(length(overlap) != 0){
      hg38_to_mm10[queryHits(overlap)]$br_score <- list[i]  
      
      match <- match(hg38_tombo_0109$site_pos_hg38, hg38_to_mm10[queryHits(overlap)]$site_pos_hg38) 
      qhits <- which(match != 'NA') 
      shits <- na.omit(match)
      
      hg38_tombo_0109$br_score[qhits] <- list[i]
      hg38_to_mm10 <- hg38_to_mm10[-queryHits(overlap)] #从hg38_to_mm10中移除匹配上的位点
    }
    print(i)
    
  }
  
  
}

#-------------------------------------------------------------------------------
saveRDS(hg38_tombo_0109,'../file/hg38_tombo_0109_br.rds')


