################################################################################
##############################      加载R包      ###############################
################################################################################
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
library(e1071)
library(caret)
library(m6ALogisticModel)
library(ggplot2) #作图
library(gridExtra) #作图
################################################################################
#############################     设置工作目录      ############################
################################################################################
setwd('/Users/huangjiaming/Desktop/conservation_demo/file')
################################################################################
#################     randomly selected As vs mouse m6A       ##################
################################################################################
m6A_SB_mm10 <- readRDS("./m6A_SB_mm10.rds") #鼠保守m6A数据库

#randomly extract As from the same transcript which m6A located #以下这部分是为了从已知甲基化位点对应的转录本中提取出所有具有m6A特征的位点，和鼠的甲基化位点取交集，看看交集分布
  #如果个数比[已知甲基化位点和鼠的甲基化位点]的交集少，则说明已知甲基化位点和鼠的甲基化的交集数目是显著的（仔细阅读results部分第一段）

hg19InMm10_m6AConservationInfo<-readRDS('./hg19InMm10_m6AConservationInfo_updated_3.rds')
#hg19InMm10_m6AConservationInfo：已知的保守位点在人当中的数据库（确认保守），看看用我们的方式是不是比别的未确定保守的多
#hg19InMm10_m6AConservationInfo：human m6A RNA methylation sites
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
bm6a <- subsetByOverlaps(transcripts(txdb),hg19InMm10_m6AConservationInfo) 

cm6a <- m6ALogisticModel::sample_sequence("DRACH",bm6a,Hsapiens,Fixed = F)
saveRDS(cm6a,'./cm6a.rds')

cm6a <- readRDS('./cm6a.rds')
dm6a <- m6ALogisticModel::sample_sequence("A",cm6a-2,Hsapiens,Fixed = T)

#cm6a数据太大，内存溢出，分段处理-----------------------------------------------
cm6a <- readRDS('./cm6a.rds')
cm6a_1 <- cm6a[1:10000000,]
saveRDS(cm6a_1,'./cm6a_1.rds')
dm6a_1 <- m6ALogisticModel::sample_sequence("A",cm6a_1-2,Hsapiens,Fixed = T)
saveRDS(dm6a_1,'./dm6a_1.rds')

cm6a <- readRDS('./cm6a.rds')
cm6a_2 <- cm6a[10000001:20000000,]
saveRDS(cm6a_2,'./cm6a_2.rds')
dm6a_2 <- m6ALogisticModel::sample_sequence("A",cm6a_2-2,Hsapiens,Fixed = T)
saveRDS(dm6a_2,'./dm6a_2.rds')

cm6a <- readRDS('./cm6a.rds')
cm6a_3 <- cm6a[20000001:22238032,]
saveRDS(cm6a_3,'./cm6a_3.rds')
dm6a_3 <- m6ALogisticModel::sample_sequence("A",cm6a_3-2,Hsapiens,Fixed = T)
saveRDS(dm6a_3,'./dm6a_3.rds')

dm6a_1 <- readRDS('./dm6a_1.rds')
dm6a_2 <- readRDS('./dm6a_2.rds')
dm6a_3 <- readRDS('./dm6a_3.rds')

dm6a <- c(dm6a_1, dm6a_2, dm6a_3)
saveRDS(dm6a,'./dm6a.rds')
#-------------------------------------------------------------------------------
dm6a <- readRDS('./dm6a.rds')

dm6a$ref <- DNAStringSet(Views(Hsapiens,dm6a)) #提取满足甲基化位点条件坐标对应的碱基
table(dm6a$ref)

result <- NA
for (i in 1:1000) {
  randomA <- dm6a[sample(seq_along(dm6a),length(hg19InMm10_m6AConservationInfo))]       
  #从dm6a中随机取和已知甲基化位点数一样数量的位点
  chain <- import.chain("./hg19ToMm10.over.chain") #人hg19基因组坐标和鼠mm10坐标的转换文件
  hg19_to_mm10 <- unlist(liftOver(randomA,chain)) #将位点转为mm10的坐标
  result2 <- length(queryHits(findOverlaps(hg19_to_mm10,m6A_SB_mm10))) #和鼠的甲基化位点取交集,得到交集个数
  result <- c(result,result2) #合并所有随机后的交集结果
  print(paste0('round_',i,':',result2))
}
result <- result[-1] #第一个为NA，移除

saveRDS(result,'./RandomTestA_result.rds')

################################################################################
################################     作图      #################################
################################################################################
RandomTestA_result <- readRDS('./RandomTestA_result.rds')
result_mean <- mean(RandomTestA_result) 

#ggplot: geom_histogram只需要一个参数，即位点交集数，函数会自动计算对应位点下的count值（纵坐标）
sp <- ggplot(data=as.data.frame(RandomTestA_result),aes(RandomTestA_result))                                                  + geom_histogram(color='black',fill='white',binwidth = 10)                                                                    + scale_x_continuous(name = 'Number of human As overlapped with mouse m6As', breaks=seq(400,600,50), limits = c(400,900),expand=c(0,0))                                                                                                              + scale_y_continuous(name = 'Frequency', breaks=seq(0,250,50), limits = c(0,250),expand=c(0,0))                               + geom_line(stat='density')+theme_classic()
#控制颜色和每一个分组的大
#breaks 控制坐标的标签展示方式

a <- sp + theme( #theme: 控制图形展示的样式，例如背景色等
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

ggsave('./randomA.pdf',a,height = 4,width = 5)
#实际距离----------------------------------------------------------------------------------------------------------
sp <- ggplot(data=as.data.frame(RandomTestA_result),aes(RandomTestA_result))                                        + geom_histogram(color='black',fill='white',binwidth = 10)                                                          + scale_x_continuous(name = 'Number of human As overlapped with mouse m6As', breaks=seq(400,600,50), limits = c(400,22800),expand=c(0,0))                                                                                              + scale_y_continuous(name = 'Frequency', breaks=seq(0,250,50), limits = c(0,250),expand=c(0,0))                     + geom_line(stat='density')+theme_classic()

a <- sp + theme( #theme: 控制图形展示的样式，例如背景色等
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.grid.minor = element_blank(), 
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA))

a <- a + geom_segment(aes(x=22359, xend = 22359 , y=20, yend = 0), size=0.5,arrow = arrow(length = unit(0.2,"cm")),angle = 180,colour = "red") #但如果画这个话你要把前面x的范围改一下, sp limits (x) 限定了范围，改成比22359大的值

ggsave('./randomA.pdf',a,height = 4,width = 5)
#------------------------------------------------------------------------------------------------------------------

