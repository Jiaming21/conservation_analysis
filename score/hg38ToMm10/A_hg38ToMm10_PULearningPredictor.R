################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/score/hg38ToMm10/file_A_hg38ToMm10_PULearningPredictor')
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

match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA')
shits <- na.omit(match)
hg38_tombo_0109$status <- 'not conserved'
hg38_tombo_0109$status[qhits] <- 'conserved'

#positive sample
hg38_tombo_0109_positive <- hg38_tombo_0109[hg38_tombo_0109$status == 'conserved'] #2582行
saveRDS(hg38_tombo_0109_positive,'_pos_gr.rds')
#negative sample
hg38_tombo_0109_negative <- hg38_tombo_0109[hg38_tombo_0109$status == 'not conserved'] #151617行
saveRDS(hg38_tombo_0109_negative,'_neg_gr.rds')
