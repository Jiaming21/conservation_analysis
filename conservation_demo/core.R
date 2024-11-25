################################################################################
##############################      加载R包      ###############################
################################################################################
library(jsonlite)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
library(e1071)
library(ROCR)
library(pROC)
library(dplyr)
library(caret)
library(m6ALogisticModel)
################################################################################
#############################     设置工作目录      ############################
################################################################################
setwd('/Users/huangjiaming/Desktop/conservation_demo/file')
################################################################################
##############################      样本导入      ##############################
################################################################################
sample_input<-read.table('./sample.txt')
################################################################################
##############################      m6A提取      ###############################
################################################################################
#输入sample_input: sample_input是输入样本
source('../part1_sample_input.R') 
#输出PIANO_INPUT: 符合“DRACH”的单个碱基位置
################################################################################
#######################      MAP到已知m6A数据库      ###########################
################################################################################
#输入PIANO_INPUT: 符合“DRACH”的单个碱基位置
source('../part2_mapWithDatabase.R') 
#database_overlap: 输入数据中已发现的m6A位点
#输出newdata: 新发现的m6A
################################################################################
###############################      打分      #################################
################################################################################
#输入newdata: 新发现的m6A
source('../part3.1_scoring_liftover.R')
#输出:
#1. hg19_to_mm10: 转换到鼠的位点
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)
#3. newdata: 新发现的m6A(含metadata)

source('../part3.1.1_tissue-specific.R')

#输入:
#1. hg19_to_mm10: 转换到鼠的位点
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)
#3. newdata: 新发现的m6A(含metadata)
source('../part3.2_scoring_base-resolution.R')
#输出:
#1. hg19_to_mm10: 转换到鼠的位点(删了所有“广义保守”的行)
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)
#3. newdata: 新发现的m6A(含更多metadata)

#输入:
#1. hg19_to_mm10: 转换到鼠的位点(删了所有“广义保守”的行)
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)
#3. newdata: 新发现的m6A(含更多metadata)
source('../part3.3_scoring_tissue-specific.R')
#输出:
#1. hg19_to_mm10: 转换到鼠的位点(删了所有“广义保守”的行)(此步没改)
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新发现的m6A(含更多metadata)(增加了brts_score)

#输入:
#1. hg19_to_mm10: 转换到鼠的位点(删了所有“广义保守”的行)(此步没改)
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新发现的m6A(含更多metadata)(增加了"brts_score")
source('../part3.4_scoring_numberOfDatasetSupported.R')
#输出: 
#1. hg19_to_mm10: 转换到鼠的位点(删了所有“广义保守”的行)(此步没改)
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新发现的m6A(含更多metadata)(增加了"nodsHg19_score"和"nodsmm10_score")

#输入:
#1. hg19_to_mm10: 转换到鼠的位点(删了所有“广义保守”的行)(此步没改)
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新发现的m6A(含更多metadata)(增加了"nodsHg19_score"和"nodsmm10_score")
source('../part3.5_scoring_sequenceAnalysis.R')
#输出:
#1. hg19_to_mm10: 全新的对应碱基都为A，含有"A_11bp"，"A_mm10"，"A_11bp_mm10"，"status"，"score"
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新加"sequence_score"

#输入:
#1. hg19_to_mm10: 全新的对应碱基都为A，含有"A_11bp"，"A_mm10"，"A_11bp_mm10"，"status"，"score"
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新加"sequence_score"
source('../part3.6_scoring_PUlearningPredictor.R')
#输出:




#输入:
source('../part3.7_scoring_phastconsScore.R')
#输出:






