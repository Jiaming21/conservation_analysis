################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/score/hg38ToMm10/file_A_hg38ToMm10_numberOfDatasetSupported')
################################################################################
####################      Number of dataset supported      #####################
################################################################################
#人位点研究数打分
hg38_tombo_0109 <- readRDS('./hg38_tombo_0109.rds')
hg38_tombo_0109 <- hg38_tombo_0109[which(vcountPattern('single base',hg38_tombo_0109$resolution)>=1)]
hg38_tombo_0109 <- hg38_tombo_0109[which(hg38_tombo_0109$ref_base=='A')]
hg38_tombo_0109$NumberOfExperimentSupported_hg38 <- vcountPattern(';',hg38_tombo_0109$source) + 1
hg38_tombo_0109$nodsHg38_score <- 0

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 == 1)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 == 1)]$nodsHg38_score <- 0.2
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 == 2)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 == 2)]$nodsHg38_score <- 0.4
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 >= 3 & hg38_tombo_0109$NumberOfExperimentSupported <= 4)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 >= 3 & hg38_tombo_0109$NumberOfExperimentSupported <= 4)]$nodsHg38_score <- 0.6
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 >= 5 & hg38_tombo_0109$NumberOfExperimentSupported_hg38 <= 6)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 >= 5 & hg38_tombo_0109$NumberOfExperimentSupported_hg38 <= 6)]$nodsHg38_score <- 0.8
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 > 6)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_hg38 > 6)]$nodsHg38_score <- 1
}
#-------------------------------------------------------------------------------
#鼠位点研究数打分
hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))
chain <- import.chain('hg38ToMm10.over.chain')
hg38_to_mm10 <- unlist(liftOver(hg38_tombo_0109,chain))
hg38_to_mm10$site_pos_mm10 <- paste0(seqnames(hg38_to_mm10),':',start(hg38_to_mm10),' ',strand(hg38_to_mm10))

mm10_tombo_0109 <- readRDS('./mm10_tombo_0109.rds')
mm10_tombo_0109 <- mm10_tombo_0109[which(vcountPattern('single base',mm10_tombo_0109$resolution)>=1)]
mm10_tombo_0109 <- mm10_tombo_0109[which(mm10_tombo_0109$ref_base=='A')]
mm10_tombo_0109$NumberOfExperimentSupported_mm10 <- vcountPattern(';',mm10_tombo_0109$source) + 1

mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))

match <- match(hg38_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_to_mm10$NumberOfExperimentSupported_mm10 <- 0
hg38_to_mm10$NumberOfExperimentSupported_mm10[qhits] <- mm10_tombo_0109$NumberOfExperimentSupported_mm10[shits]

match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$NumberOfExperimentSupported_mm10 <- 0
hg38_tombo_0109$NumberOfExperimentSupported_mm10[qhits] <- hg38_to_mm10$NumberOfExperimentSupported_mm10[shits]


hg38_tombo_0109$nodsMm10_score <- 0
if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 == 1)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 == 1)]$nodsMm10_score <- 0.2
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 == 2)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 == 2)]$nodsMm10_score <- 0.4
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 >= 3 & hg38_tombo_0109$NumberOfExperimentSupported_mm10 <= 4)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 >= 3 & hg38_tombo_0109$NumberOfExperimentSupported_mm10 <= 4)]$nodsMm10_score <- 0.6
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 >= 5 & hg38_tombo_0109$NumberOfExperimentSupported_mm10 <= 6)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 >= 5 & hg38_tombo_0109$NumberOfExperimentSupported_mm10 <= 6)]$nodsMm10_score <- 0.8
}

if(length(which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 > 6)) != 0){
  hg38_tombo_0109[which(hg38_tombo_0109$NumberOfExperimentSupported_mm10 > 6)]$nodsMm10_score <- 1
}
#-------------------------------------------------------------------------------
saveRDS(hg38_tombo_0109,'../file/hg38_tombo_0109_nodp.rds')


