#! /usr/bin/env Rscript
################################################################################
##############################      加载R包      ###############################
################################################################################
library(GenomicRanges)

library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(m6ALogisticModel)
library(ggplot2) #作图
library(gridExtra) #作图
################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/validation/file_hg38ToMm10')
################################################################################
############################      统计显著性分析    ############################
################################################################################
#三代测序：hg38_tombo_0109 mm10_tombo_0109 #泛A修饰
#hg38_tombo_0109 <- readRDS('hg38_tombo_0109.rds')
#hg38_tombo_0109 <- hg38_tombo_0109[which(vcountPattern('single base',hg38_tombo_0109$resolution)>=1)]
#saveRDS(hg38_tombo_0109,'hg38_tombo_0109.rds')
#mm10_tombo_0109 <- readRDS('mm10_tombo_0109.rds')
#mm10_tombo_0109 <- mm10_tombo_0109[which(vcountPattern('single base',mm10_tombo_0109$resolution)>=1)]
#saveRDS(mm10_tombo_0109,'mm10_tombo_0109.rds')

#读取数据
hg38_tombo_0109 <- readRDS('hg38_tombo_0109.rds')
mm10_tombo_0109 <- readRDS('mm10_tombo_0109.rds')
#研究对象定为A修饰
hg38_tombo_0109 <- hg38_tombo_0109[hg38_tombo_0109$ref_base == 'A'] #154199
hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))
chain <- import.chain('hg38ToMm10.over.chain')
hg38_to_mm10 <- unlist(liftOver(hg38_tombo_0109,chain))
match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$liftOverToMm10 <- 'No'
hg38_tombo_0109$liftOverToMm10[qhits] <- 'Yes'

hg38_to_mm10$site_pos_mm10 <- paste0(seqnames(hg38_to_mm10),':',start(hg38_to_mm10),' ',strand(hg38_to_mm10))
mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))
match <- match(hg38_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_to_mm10$modification_mm10 <- '-'
hg38_to_mm10$modification_mm10[qhits] <- mm10_tombo_0109$ref_base[shits]

length(hg38_to_mm10[hg38_to_mm10$modification_mm10 == 'A']) #qhits=shits: 能成功配上的行数 #2582

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
  hg38_to_mm10 <- unlist(liftOver(randomA,chain)) #新的hg38_to_mm10（由randomA转化得来）
  match <- match(randomA$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  randomA$liftOverToMm10 <- 'No'
  randomA$liftOverToMm10[qhits] <- 'Yes'
  
  hg38_to_mm10$site_pos_mm10 <- paste0(seqnames(hg38_to_mm10),':',start(hg38_to_mm10),' ',strand(hg38_to_mm10))
  mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))
  match <- match(hg38_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  hg38_to_mm10$modification_mm10 <- '-'
  hg38_to_mm10$modification_mm10[qhits] <- mm10_tombo_0109$ref_base[shits]
  result2 <- length(hg38_to_mm10[hg38_to_mm10$modification_mm10 == 'A']) 
  result <- c(result,result2) 
  print(paste0('round_',i,':',result2))
}
result <- result[-1]
saveRDS(result,'./RandomTestA_result.rds')

################################################################################
################################     作图      #################################
################################################################################
RandomTestA_result <- readRDS('./RandomTestA_result.rds')
result_mean <- mean(RandomTestA_result) #8.198

#ggplot: geom_histogram只需要一个参数，即位点交集数，函数会自动计算对应位点下的count值（纵坐标）
sp <- ggplot(data=as.data.frame(RandomTestA_result),aes(RandomTestA_result)) + geom_histogram(color='black',fill='white',binwidth = 1) + scale_x_continuous(name = 'Number of human As overlapped with mouse As', breaks=seq(0,35,5), limits = c(-0.5,100),expand=c(0,0)) + scale_y_continuous(name = 'Frequency', breaks=seq(0,300,50), limits = c(0,300),expand=c(0,0)) + geom_line(stat='density') + theme_classic()

#控制颜色和每一个分组的大
#breaks 控制坐标的标签展示方式

a <- sp + theme( #theme: 控制图形展示的样式，例如背景色等
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

ggsave('./randomA_hg38ToMm10.pdf',a,height = 4,width = 5)


