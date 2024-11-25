#前情回顾: 
#1. hg19_to_mm10: 全新的对应碱基都为A，含有"A_11bp"，"A_mm10"，"A_11bp_mm10"，"status"，"score"
#2. m6A_SB_mm10: 鼠保守m6A数据库(一直没变)(此步没改)
#3. newdata: 新加"sequence_score"

################################################################################
######################      PU learning predictor      #########################
################################################################################
#Genomic feature generation
#hg19  
#主要是根据m6ALogisticModel的函数功能进行一些操作, 详见: https://github.com/ZW-xjtlu/m6ALogisticModel
GFgenreation_m6A <- function(data){
  analysis_data <- data
  matureSE <- SummarizedExperiment() 
  #定义一个SummarizedExperiment（SE）对象，SE是矩阵样容器，行是特征(genes, transcripts, exons, etc.)，列是样本。
  #每个对象储存着一个或多个样本的观测结果，以及用于描述特征信息(features)和样本信息 (phenotypes) 的额外元数据（metadata）
  rowRanges(matureSE) <- analysis_data  #rowRanges存储特征信息
  Additional_features_hg19 = list( #自定义特征，各种数据库的信息
    HNRNPC_eCLIP = eCLIP_HNRNPC_gr, 
    YTHDC1_TREW = YTHDC1_TREW_gr,
    YTHDF1_TREW = YTHDF1_TREW_gr,
    YTHDF2_TREW = YTHDF2_TREW_gr,
    miR_targeted_genes = miR_targeted_genes_grl,
    TargetScan = TargetScan_hg19_gr,
    Verified_miRtargets = verified_targets_gr,
    METTL3_TREW = METTL3_TREW,
    METTL14_TREW = METTL14_TREW,
    WTAP_TREW = WTAP_TREW,
    METTL16_CLIP = METTL16_CLIP,
    ALKBH5_PARCLIP = ALKBH5_PARCLIP,
    FTO_CLIP = FTO_CLIP,
    FTO_eCLIP = FTO_eCLIP
  )
  #融合各个数据库信息将各个位点相应的特征提取出来
  data_standardized <- predictors_annot(se = matureSE,#SE格式的数据
                                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, #人类基因组注释,位置信息
                                        bsgnm = Hsapiens, #人类的全基因组序列
                                        fc = fitCons.UCSC.hg19, #standardized Fitness consequences scores
                                        pc = phastCons100way.UCSC.hg19, #保守性得分,UCSC phastCons conservation scores
                                        struct_hybridize = Struc_hg19, #Hybrid zones are locations where hybrids between species, subspecies, or races are found
                                        feature_lst = Additional_features_hg19,#自定义的其他特征
                                        hk_genes_list = HK_hg19_eids, #house keeping genes list
                                        motif = c("DRACH"), #A character vector indicating the motifs centered by the modification nucleotite
                                        motif_clustering = "A", #A character vector indicating the motif used to generate the features for the clustering indexes
                                        isoform_ambiguity_method = "longest_tx",# keeps only the longest transcript as the transcript annotation
                                        genes_ambiguity_method = "average",#use the average feature entries for mapping of multiple genes.
                                        annot_clustering = matureSE,#A GRanges object to generate clustering features
                                        standardization = T) #A logical indicating whether to standardize the continous features
  
  GF <- mcols(data_standardized)
  return(GF)
}

background <- readRDS('./hg19InMm10_m6AConservationInfo_updated_3.rds')[1:2400]
newdata <- c(background,newdata) #将background和newdata合并，对结果的标准化有影响
genome(newdata) <- NA #genome 属性为空
newdata <- newdata[,1:26] #取前26列
GF_positiveAll <- GFgenreation_m6A(newdata) #取特征信息
GF_positiveAll <- as.data.frame(GF_positiveAll)

for (i in 1:ncol(GF_positiveAll)) {
  if(anyNA(GF_positiveAll[,i]) == TRUE){ #有NA的列都赋值为0
    GF_positiveAll[,i] <- 0
    print(i)
  }
}

GF_positiveAll <- GF_positiveAll[,-30] #移除motif DRACH, 都为1

#------------------------------------------------------------------------------#
#Sequence feature generation
SeqFgeneration <- function(data,GF){ #GF无用？
  UntestVarSeq <- data
  source('../method/class1.R')
  source('../method/class2.R')
  source('../method/class3.R')
  source('../method/class4.R')
  source('../method/class5.R')
  source('../method/class6.R')
  CP <- ChemicalProperty(UntestVarSeq)
  return(CP)
}
SeqFgeneration2 <- function(data,GF){ #GF无用？
  UntestVarSeq <- data
  source('../method/class1.R')
  source('../method/class2.R')
  source('../method/class3.R')
  source('../method/class4.R')
  source('../method/class5.R')
  source('../method/class6.R')
  NF <- sequenceFeatures(UntestVarSeq,NTYPE="DNA")
  return(NF)
}
sequenceFeature1 <- SeqFgeneration(as.character(newdata$reference_sequence),GF)#？GF_positiveAll #将序列中相应的碱基（ATCGNT）转换为01编码，每个碱基对应3个数字编码；第21个碱基都是A，第22个碱基都是C所以61:66编码是一样的
sequenceFeature2 <- SeqFgeneration2(as.character(newdata$reference_sequence),GF) #返回每个碱基之前有多少比例的相同碱基矩阵;第一列都是1
sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)],sequenceFeature2[,-1]) 

#------------------------------------------------------------------------------#
#特征：
FeatureBoth <- cbind(GF_positiveAll,sequenceFeature) #将所有特征信息合并
FeatureBoth <- FeatureBoth[2401:length(newdata),] #只取输入文件中特有的行
#------------------------------------------------------------------------------#
newdata <- newdata[2401:length(newdata)]
rownames(FeatureBoth) <- seq_along(newdata) #seq_along创建一个从1开始、长度为其输入值的序列,注意和seq_len区别
#------------------------------------------------------------------------------#
#best model - Caret_conservationHg19_model_9.rds
conservationHg19_model_9 <- readRDS('./Caret_conservationHg19_model_9.rds') #调已经训练好的模型
BIOtest_pred <- predict(conservationHg19_model_9, newdata = FeatureBoth, type = "prob") #输出预测概率

indxMiss <- setdiff(1:length(newdata),rownames(BIOtest_pred)) #求向量x与向量y中不同的元素(只取x中不同的元素) ,考虑可能有行无法预测，缺失
if(length(indxMiss) != 0){
  newdata <- newdata[-indxMiss,]
}

newdata$PUlearningHg19_score <- BIOtest_pred$positive #取预测positive概率为分值


#mm10
#best model - Caret_conservationMm10_model_5.rds
indx_A <- which(hg19_to_mm10$A_mm10 == 'A') #对应鼠中的位点必须是A

if(length(indx_A) != 0){
  conservationMm10_model_5 <- readRDS('./Caret_conservationMm10_model_5.rds') #鼠中训练的模型
  GFgenreation_mm10m6A <- function(data){
    analysis_data <- data
    matureSE <- SummarizedExperiment()
    rowRanges(matureSE) <- analysis_data 
    data_standardized <- predictors_annot(se = matureSE,
                                          txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, #对应的基因组换成鼠的
                                          bsgnm = Mmusculus,
                                          struct_hybridize = Struc_mm10,
                                          motif_clustering = "A",
                                          isoform_ambiguity_method = "longest_tx",
                                          genes_ambiguity_method = "average",
                                          annot_clustering = matureSE,
                                          standardization = T) 
    GF <- mcols(data_standardized)
    return(GF) 
  }
  
  background_mm10 <- readRDS("./m6A_SB_mm10.rds")[1:2400] #鼠的background甲基化位点
  mcols(background_mm10) <- mcols(hg19_to_mm10[1]) #将background_mm10的metadata column复制成hg19_to_mm10第一行，这些列暂时没有用，主要是为了使格式统一
  hg19_to_mm10 <- c(background_mm10,hg19_to_mm10) #合并行
  GF_mm10 <- GFgenreation_mm10m6A(hg19_to_mm10) #同上，得到基因组位置等特征信息
  GF_mm10 <- as.data.frame(GF_mm10)
  GF_positiveAll_mm10 <- readRDS('./GF_positiveAll_mm10.rds')[2,] #获取列的信息,数据框,可能鼠模型用到的列？
  match <- match(names(GF_mm10),names(GF_positiveAll_mm10)) #names = colnames
  GF_mm10 <- GF_mm10[,which(match != 'NA')] #去除不在GF_positiveAll_mm10.rds中的列
  
  for (i in 1:ncol(GF_mm10)) {
    if(anyNA(GF_mm10[,i]) == TRUE){
      GF_mm10[,i] <- 0
      print(i)
    }
  }
  
  sequenceFeature1 <- SeqFgeneration(as.character(DNAStringSet(Views(Mmusculus,hg19_to_mm10 + 20))),GF) #同上，取对应鼠的序列
  sequenceFeature2 <- SeqFgeneration2(as.character(DNAStringSet(Views(Mmusculus,hg19_to_mm10 + 20))),GF)
  sequenceFeature <- cbind(sequenceFeature1[,-c(61:66)],sequenceFeature2[,-1])
  featureMm10 <- cbind(GF_mm10,sequenceFeature) #合并所有特征
  featureMm10 <- featureMm10[2401:length(hg19_to_mm10),]
  hg19_to_mm10 <- hg19_to_mm10[2401:length(hg19_to_mm10)]
  rownames(featureMm10) <- seq_along(hg19_to_mm10)
  BIOtest_pred_mm10 <- predict(conservationMm10_model_5, newdata = featureMm10, type = "prob") #预测
  
  indxMiss <- setdiff(1:length(hg19_to_mm10),rownames(BIOtest_pred_mm10))
  if(length(indxMiss) != 0){
    hg19_to_mm10 <- hg19_to_mm10[-indxMiss,]
  }
  
  hg19_to_mm10$PUlearningMm10_score <- BIOtest_pred_mm10$positive
  newdata$PUlearningMm10_score <- 0
  
  match <- match(newdata$ID,hg19_to_mm10$ID)
  qhits <- which(match != 'NA')
  shits <- na.omit(match)
  hg19_to_mm10$A_11bp_mm10 <- DNAStringSet(Views( Mmusculus,hg19_to_mm10 + 5 )) #鼠的11bp的序列
  newdata$PUlearningMm10_score[qhits] <- hg19_to_mm10$PUlearningMm10_score[shits] #对应位置赋值
  newdata$A11bp_hmm10 <- '-'
  newdata$base_mm10 <- 'other'
  newdata$A11bp_hmm10[qhits] <- hg19_to_mm10$A_11bp_mm10[shits]
  newdata$base_mm10[qhits] <- 'A' #补充列的信息
}else{
  newdata$PUlearningMm10_score <- 0
  newdata$A11bp_hmm10 <- '-'
  newdata$base_mm10 <- 'other'
}


indx <- which(newdata$conserved_status == 'conserved')
if(length(indx) != 0){
  newdata[indx]$PUlearningHg19_score <- 1 #保守位置的分值直接设置为1
  newdata[indx]$PUlearningMm10_score <- 1
}





