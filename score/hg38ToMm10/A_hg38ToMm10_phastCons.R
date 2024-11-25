################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/score/hg38ToMm10/file_A_hg38ToMm10_phastCons')
################################################################################
##########################      PhastCons score      ###########################
################################################################################
library(phastCons100way.UCSC.hg38)
gsco <- phastCons100way.UCSC.hg38
hg38_tombo_0109 <- readRDS('./hg38_tombo_0109.rds')
hg38_tombo_0109 <- hg38_tombo_0109[which(vcountPattern('single base',hg38_tombo_0109$resolution)>=1)]
hg38_tombo_0109 <- hg38_tombo_0109[which(hg38_tombo_0109$ref_base=='A')]
hg38_tombo_0109 <- gscores(gsco,hg38_tombo_0109,pop="DP2")
indx_na <- as.numeric(attr(na.omit(hg38_tombo_0109$DP2),"na.action"))
hg38_tombo_0109$DP2[indx_na] <- 0