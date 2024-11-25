#! /usr/bin/env Rscript
################################################################################
##############################      加载R包      ###############################
################################################################################
library(GenomicRanges)

library(BSgenome.Sscrofa.UCSC.susScr11)
library(TxDb.Sscrofa.UCSC.susScr11.refGene)

library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

library(m6ALogisticModel)
library(ggplot2) #作图
library(gridExtra) #作图
################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/validation/file_sus11ToMm10')
################################################################################
############################      统计显著性分析    ############################
################################################################################
#读取数据
#sus11 <- read.table('SusScrofaunknown_site.txt',row.names = NULL)
#colnames(sus11) <- colnames(sus11)[-1]
#sus11[which(sus11$seqnames == 'NC_000845.1'),]$seqnames <- 'chrM'
#sus11[which(sus11$seqnames == 'NC_010443.5'),]$seqnames <- 'chr1'
#sus11[which(sus11$seqnames == 'NC_010444.4'),]$seqnames <- 'chr2'
#sus11[which(sus11$seqnames == 'NC_010445.4'),]$seqnames <- 'chr3'
#sus11[which(sus11$seqnames == 'NC_010446.5'),]$seqnames <- 'chr4'
#sus11[which(sus11$seqnames == 'NC_010447.5'),]$seqnames <- 'chr5'
#sus11[which(sus11$seqnames == 'NC_010448.4'),]$seqnames <- 'chr6'
#sus11[which(sus11$seqnames == 'NC_010449.5'),]$seqnames <- 'chr7'
#sus11[which(sus11$seqnames == 'NC_010450.4'),]$seqnames <- 'chr8'
#sus11[which(sus11$seqnames == 'NC_010451.4'),]$seqnames <- 'chr9'
#sus11[which(sus11$seqnames == 'NC_010452.4'),]$seqnames <- 'chr10'
#sus11[which(sus11$seqnames == 'NC_010453.5'),]$seqnames <- 'chr11'
#sus11[which(sus11$seqnames == 'NC_010454.4'),]$seqnames <- 'chr12'
#sus11[which(sus11$seqnames == 'NC_010455.5'),]$seqnames <- 'chr13'
#sus11[which(sus11$seqnames == 'NC_010456.5'),]$seqnames <- 'chr14'
#sus11[which(sus11$seqnames == 'NC_010457.5'),]$seqnames <- 'chr15'
#sus11[which(sus11$seqnames == 'NC_010458.4'),]$seqnames <- 'chr15'
#sus11[which(sus11$seqnames == 'NC_010459.5'),]$seqnames <- 'chr17'
#sus11[which(sus11$seqnames == 'NC_010460.4'),]$seqnames<- 'chr18'
#sus11[which(sus11$seqnames == 'NC_010461.5'),]$seqnames <- 'chrX'
#sus11[which(sus11$seqnames == 'NC_010462.3'),]$seqnames <- 'chrY'
#sus11[which(sus11$seqnames == 'NW_018084833.1'),]$seqnames <- 'chrUn_NW_018084833v1'
#sus11[which(sus11$seqnames == 'NW_018084901.1'),]$seqnames <- 'chrUn_NW_018084901v1'
#sus11[which(sus11$seqnames == 'NW_018084979.1'),]$seqnames <- 'chrUn_NW_018084979v1'
#sus11[which(sus11$seqnames == 'NW_018085108.1'),]$seqnames <- 'chrUn_NW_018085108v1'
#sus11[which(sus11$seqnames == 'NW_018085356.1'),]$seqnames <- 'chrUn_NW_018085356v1'
#ref_base <- sus11$ref_base
#sus11 <- GRanges(seqnames = sus11$seqnames,
#                 ranges = IRanges(start = sus11$start,
#                                  end = sus11$end,
#                                  width = 1),
#                 strand = sus11$strand)
#sus11$ref_base <- ref_base
#saveRDS(sus11,'SusScrofaunknown_site.rds')

sus11 <- readRDS('SusScrofaunknown_site.rds')
mm10_tombo_0109 <- readRDS('mm10_tombo_0109.rds')

#研究对象定为A修饰
sus11 <- sus11[sus11$ref_base == 'A'] #23060
sus11$site_pos_sus11 <- paste0(seqnames(sus11),':',start(sus11),' ',strand(sus11))
chain <- import.chain('susScr11ToMm10.over.chain')
sus11_to_mm10 <- unlist(liftOver(sus11,chain))
match <- match(sus11$site_pos_sus11,sus11$site_pos_sus11)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
sus11$liftOverToMm10 <- 'No'
sus11$liftOverToMm10[qhits] <- 'Yes'

sus11_to_mm10$site_pos_mm10 <- paste0(seqnames(sus11_to_mm10),':',start(sus11_to_mm10),' ',strand(sus11_to_mm10))
mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))
match <- match(sus11_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
sus11_to_mm10$modification_mm10 <- '-'
sus11_to_mm10$modification_mm10[qhits] <- mm10_tombo_0109$ref_base[shits]

length(sus11_to_mm10[sus11_to_mm10$modification_mm10 == 'A']) #qhits=shits: 能成功配上的行数 #40

#-----------------------------------------------------------------------------------------------
txdb <- TxDb.Sscrofa.UCSC.susScr11.refGene

deleted <- subsetByOverlaps(transcripts(txdb),sus11) 
A <- m6ALogisticModel::sample_sequence("A",deleted,Sscrofa,Fixed = T)
saveRDS(A,'./A.rds')
A <- readRDS('./A.rds')

result <- NA
for (i in 1:1000) {
  randomA <- A[sample(seq_along(A),length(sus11))]
  randomA$site_pos_sus11 <- paste0(seqnames(randomA),':',start(randomA),' ',strand(randomA))
  sus11_to_mm10 <- unlist(liftOver(randomA,chain)) #新的hg38_to_mm10（由randomA转化得来）
  match <- match(randomA$site_pos_sus11,sus11_to_mm10$site_pos_sus11)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  randomA$liftOverToMm10 <- 'No'
  randomA$liftOverToMm10[qhits] <- 'Yes'
  
  sus11_to_mm10$site_pos_mm10 <- paste0(seqnames(sus11_to_mm10),':',start(sus11_to_mm10),' ',strand(sus11_to_mm10))
  mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))
  match <- match(sus11_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
  qhits <- which(match != 'NA') 
  shits <- na.omit(match)
  sus11_to_mm10$modification_mm10 <- '-'
  sus11_to_mm10$modification_mm10[qhits] <- mm10_tombo_0109$ref_base[shits]
  result2 <- length(sus11_to_mm10[sus11_to_mm10$modification_mm10 == 'A']) 
  result <- c(result,result2) 
  print(paste0('round_',i,':',result2))
}
result <- result[-1]
saveRDS(result,'./RandomTestA_result.rds')

################################################################################
################################     作图      #################################
################################################################################
RandomTestA_result <- readRDS('./RandomTestA_result.rds')
result_mean <- mean(RandomTestA_result) #11

#ggplot: geom_histogram只需要一个参数，即位点交集数，函数会自动计算对应位点下的count值（纵坐标）
sp <- ggplot(data=as.data.frame(RandomTestA_result),aes(RandomTestA_result)) + geom_histogram(color='black',fill='white',binwidth = 1) + scale_x_continuous(name = 'Number of pig As overlapped with mouse As', breaks=seq(0,35,5), limits = c(-0.5,100),expand=c(0,0)) + scale_y_continuous(name = 'Frequency', breaks=seq(0,300,50), limits = c(0,300),expand=c(0,0)) + geom_line(stat='density') + theme_classic()

#控制颜色和每一个分组的大
#breaks 控制坐标的标签展示方式

a <- sp + theme( #theme: 控制图形展示的样式，例如背景色等
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

ggsave('./randomA_sus11ToMm10.pdf',a,height = 4,width = 5)
