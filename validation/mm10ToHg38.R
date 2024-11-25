#! /usr/bin/env Rscript
################################################################################
##############################      加载R包      ###############################
################################################################################
library(GenomicRanges)

library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(m6ALogisticModel)
library(ggplot2) #作图
library(gridExtra) #作图
################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/validation/file_mm10ToHg38')
################################################################################
############################      统计显著性分析    ############################
################################################################################
#三代测序：hg38_tombo_0109 mm10_tombo_0109 #泛A修饰
#读取数据
hg38_tombo_0109 <- readRDS('hg38_tombo_0109.rds')
mm10_tombo_0109 <- readRDS('mm10_tombo_0109.rds')
#研究对象定为A修饰
mm10_tombo_0109 <- mm10_tombo_0109[mm10_tombo_0109$ref_base == 'A'] #19633
mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))
chain <- import.chain('mm10ToHg38.over.chain')
mm10_to_hg38 <- unlist(liftOver(mm10_tombo_0109,chain))
match <- match(mm10_tombo_0109$site_pos_mm10,mm10_to_hg38$site_pos_mm10)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
mm10_tombo_0109$liftOverToHg38 <- 'No'
mm10_tombo_0109$liftOverToHg38[qhits] <- 'Yes'

mm10_to_hg38$site_pos_hg38 <- paste0(seqnames(mm10_to_hg38),':',start(mm10_to_hg38),' ',strand(mm10_to_hg38))
hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))
match <- match(mm10_to_hg38$site_pos_hg38,hg38_tombo_0109$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
mm10_to_hg38$modification_hg38 <- '-'
mm10_to_hg38$modification_hg38[qhits] <- hg38_tombo_0109$ref_base[shits]

length(mm10_to_hg38[mm10_to_hg38$modification_hg38 == 'A']) #qhits=shits: 能成功配上的行数 #2565

#-----------------------------------------------------------------------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

deleted <- subsetByOverlaps(transcripts(txdb),mm10_tombo_0109) 
A <- m6ALogisticModel::sample_sequence("A",deleted,Mmusculus,Fixed = T)
saveRDS(A,'./A.rds')
A <- readRDS('./A.rds')

result <- NA
for (i in 1:1000) {
  randomA <- A[sample(seq_along(A),length(mm10_tombo_0109))]
  randomA$site_pos_mm10 <- paste0(seqnames(randomA),':',start(randomA),' ',strand(randomA))
  mm10_to_hg38 <- unlist(liftOver(randomA,chain)) #新的mm10_to_hg38（由randomA转化得来）
  match <- match(randomA$site_pos_mm10,mm10_to_hg38$site_pos_mm10)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  randomA$liftOverToHg38 <- 'No'
  randomA$liftOverToHg38[qhits] <- 'Yes'
  
  mm10_to_hg38$site_pos_hg38 <- paste0(seqnames(mm10_to_hg38),':',start(mm10_to_hg38),' ',strand(mm10_to_hg38))
  hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))
  match <- match(mm10_to_hg38$site_pos_hg38,hg38_tombo_0109$site_pos_hg38)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  mm10_to_hg38$modification_hg38 <- '-'
  mm10_to_hg38$modification_hg38[qhits] <- hg38_tombo_0109$ref_base[shits]
  result2 <- length(mm10_to_hg38[mm10_to_hg38$modification_hg38 == 'A']) 
  result <- c(result,result2) 
  print(paste0('round_',i,':',result2))
}
result <- result[-1]
saveRDS(result,'./RandomTestA_result.rds')

################################################################################
################################     作图      #################################
################################################################################
RandomTestA_result <- readRDS('./RandomTestA_result.rds')
result_mean <- mean(RandomTestA_result) #11.001

#ggplot: geom_histogram只需要一个参数，即位点交集数，函数会自动计算对应位点下的count值（纵坐标）
sp <- ggplot(data=as.data.frame(RandomTestA_result),aes(RandomTestA_result)) + geom_histogram(color='black',fill='white',binwidth = 1) + scale_x_continuous(name = 'Number of mouse As overlapped with human As', breaks=seq(0,35,5), limits = c(-0.5,100),expand=c(0,0)) + scale_y_continuous(name = 'Frequency', breaks=seq(0,300,50), limits = c(0,300),expand=c(0,0)) + geom_line(stat='density') + theme_classic()

#控制颜色和每一个分组的大
#breaks 控制坐标的标签展示方式

a <- sp + theme( #theme: 控制图形展示的样式，例如背景色等
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

ggsave('./randomA_mm10ToHg38.pdf',a,height = 4,width = 5)


