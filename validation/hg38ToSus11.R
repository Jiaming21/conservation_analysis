#! /usr/bin/env Rscript
################################################################################
##############################      加载R包      ###############################
################################################################################
library(GenomicRanges)

library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(BSgenome.Sscrofa.UCSC.susScr11)
library(TxDb.Sscrofa.UCSC.susScr11.refGene)

library(m6ALogisticModel)
library(ggplot2) #作图
library(gridExtra) #作图
################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/validation/file_hg38ToSus11')
################################################################################
############################      统计显著性分析    ############################
################################################################################
#读取数据
hg38_tombo_0109 <- readRDS('hg38_tombo_0109.rds')
sus11 <- readRDS('SusScrofaunknown_site.rds')
#研究对象定为A修饰
hg38_tombo_0109 <- hg38_tombo_0109[hg38_tombo_0109$ref_base == 'A'] #154199
hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))
chain <- import.chain('hg38ToSusScr11.over.chain')
hg38_to_sus11 <- unlist(liftOver(hg38_tombo_0109,chain))
match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_sus11$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$liftOverToSus11 <- 'No'
hg38_tombo_0109$liftOverToSus11[qhits] <- 'Yes'

hg38_to_sus11$site_pos_sus11 <- paste0(seqnames(hg38_to_sus11),':',start(hg38_to_sus11),' ',strand(hg38_to_sus11))
sus11$site_pos_sus11 <- paste0(seqnames(sus11),':',start(sus11),' ',strand(sus11))
match <- match(hg38_to_sus11$site_pos_sus11,sus11$site_pos_sus11)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_to_sus11$modification_sus11 <- '-'
hg38_to_sus11$modification_sus11[qhits] <- sus11$ref_base[shits]

length(hg38_to_sus11[hg38_to_sus11$modification_sus11 == 'A']) #qhits=shits: 能成功配上的行数 #400

#-----------------------------------------------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

deleted <- subsetByOverlaps(transcripts(txdb),hg38_tombo_0109) 
A <- m6ALogisticModel::sample_sequence("A",deleted,Hsapiens,Fixed = T)
saveRDS(A,'./A.rds')
A <- readRDS('./A.rds')

result <- NA
for (i in 1:1000) {
  randomA <- A[sample(seq_along(A),length(hg38_tombo_0109))]
  randomA$site_pos_hg38 <- paste0(seqnames(randomA),':',start(randomA),' ',strand(randomA))
  hg38_to_sus11 <- unlist(liftOver(randomA,chain)) #新的hg38_to_mm10（由randomA转化得来）
  match <- match(randomA$site_pos_hg38,hg38_to_sus11$site_pos_hg38)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  randomA$liftOverToSus11 <- 'No'
  randomA$liftOverToSus11[qhits] <- 'Yes'
  
  hg38_to_sus11$site_pos_sus11 <- paste0(seqnames(hg38_to_sus11),':',start(hg38_to_sus11),' ',strand(hg38_to_sus11))
  sus11$site_pos_sus11 <- paste0(seqnames(sus11),':',start(sus11),' ',strand(sus11))
  match <- match(hg38_to_sus11$site_pos_sus11,sus11$site_pos_sus11)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  hg38_to_sus11$modification_sus11 <- '-'
  hg38_to_sus11$modification_sus11[qhits] <- sus11$ref_base[shits]
  result2 <- length(hg38_to_sus11[hg38_to_sus11$modification_sus11 == 'A']) 
  result <- c(result,result2) 
  print(paste0('round_',i,':',result2))
}
result <- result[-1]
saveRDS(result,'./RandomTestA_result.rds')

################################################################################
################################     作图      #################################
################################################################################
RandomTestA_result <- readRDS('./RandomTestA_result.rds')
result_mean <- mean(RandomTestA_result) #14.542

#ggplot: geom_histogram只需要一个参数，即位点交集数，函数会自动计算对应位点下的count值（纵坐标）
sp <- ggplot(data=as.data.frame(RandomTestA_result),aes(RandomTestA_result)) + geom_histogram(color='black',fill='white',binwidth = 1) + scale_x_continuous(name = 'Number of human As overlapped with pig As', breaks=seq(0,35,5), limits = c(-0.5,100),expand=c(0,0)) + scale_y_continuous(name = 'Frequency', breaks=seq(0,300,50), limits = c(0,300),expand=c(0,0)) + geom_line(stat='density') + theme_classic()

#控制颜色和每一个分组的大
#breaks 控制坐标的标签展示方式

a <- sp + theme( #theme: 控制图形展示的样式，例如背景色等
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

ggsave('./randomA_hg38ToSus11.pdf',a,height = 4,width = 5)

