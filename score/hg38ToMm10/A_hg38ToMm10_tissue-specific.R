################################################################################
#############################      设置工作目录      ###########################
################################################################################
setwd('/home/jiaming/score/hg38ToMm10/file_A_hg38ToMm10_tissue-specific')
################################################################################
##########################      Tissue specific      ###########################
################################################################################
#人细胞信息从二代转移到三代
hg38_tombo_0109 <- readRDS('./hg38_tombo_0109.rds')
hg38_tombo_0109 <- hg38_tombo_0109[which(vcountPattern('single base',hg38_tombo_0109$resolution)>=1)]
hg38_tombo_0109 <- hg38_tombo_0109[which(hg38_tombo_0109$ref_base=='A')]
hg38_tombo_0109$site_pos_hg38 <- paste0(seqnames(hg38_tombo_0109),':',start(hg38_tombo_0109),' ',strand(hg38_tombo_0109))

hg38_ngsinfo <- readRDS('./hg38_ngsinfo.rds')
hg38_ngsinfo <- hg38_ngsinfo[which(vcountPattern('single base',hg38_ngsinfo$resolution)>=1)]
hg38_ngsinfo$site_pos_hg38 <- paste0(seqnames(hg38_ngsinfo),':',start(hg38_ngsinfo),' ',strand(hg38_ngsinfo))

match <- match(hg38_tombo_0109$site_pos_hg38,hg38_ngsinfo$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$cell_hg38[qhits] <- hg38_ngsinfo$cell[shits] 

#鼠细胞信息从二代转移到三代
mm10_tombo_0109 <- readRDS('./mm10_tombo_0109.rds')
mm10_tombo_0109 <- mm10_tombo_0109[which(vcountPattern('single base',mm10_tombo_0109$resolution)>=1)]
mm10_tombo_0109 <- mm10_tombo_0109[which(mm10_tombo_0109$ref_base=='A')]
mm10_tombo_0109$site_pos_mm10 <- paste0(seqnames(mm10_tombo_0109),':',start(mm10_tombo_0109),' ',strand(mm10_tombo_0109))

mm10_ngsinfo <- readRDS('./mm10_ngsinfo.rds')
mm10_ngsinfo <- mm10_ngsinfo[which(vcountPattern('single base',mm10_ngsinfo$resolution)>=1)]
mm10_ngsinfo$site_pos_mm10 <- paste0(seqnames(mm10_ngsinfo),':',start(mm10_ngsinfo),' ',strand(mm10_ngsinfo))

match <- match(mm10_tombo_0109$site_pos_mm10,mm10_ngsinfo$site_pos_mm10)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
mm10_tombo_0109$cell_mm10[qhits] <- mm10_ngsinfo$cell[shits] 

#-------------------------------------------------------------------------------
#是否可以转换到鼠坐标
chain <- import.chain('hg38ToMm10.over.chain')
hg38_to_mm10 <- unlist(liftOver(hg38_tombo_0109,chain))
match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$liftOverToMm10 <- 'No'
hg38_tombo_0109$liftOverToMm10[qhits] <- 'Yes'

#信息从鼠转移到鼠(人)
hg38_to_mm10$site_pos_mm10 <- paste0(seqnames(hg38_to_mm10),':',start(hg38_to_mm10),' ',strand(hg38_to_mm10))
match <- match(hg38_to_mm10$site_pos_mm10,mm10_tombo_0109$site_pos_mm10)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_to_mm10$cell_mm10 <- '-'
hg38_to_mm10$cell_mm10[qhits] <- mm10_tombo_0109$cell_mm10[shits]

#信息从鼠转移到人
match <- match(hg38_tombo_0109$site_pos_hg38,hg38_to_mm10$site_pos_hg38)
qhits <- which(match != 'NA') 
shits <- na.omit(match)
hg38_tombo_0109$cell_mm10 <- '-'
hg38_tombo_0109$cell_mm10[qhits] <- hg38_to_mm10$cell_mm10[shits]


hg38_tombo_0109$brts_score <- 0
list <- c("brain","kidney","liver","ESC")
for (i in 1:4) {
  list_hg38 <- grep(list[i],hg38_tombo_0109$cell_hg38) #grep：查找目标字符串或字符串向量中是否包含目标子串,返回位置信息index
  list_mm10 <- grep(list[i],hg38_tombo_0109$cell_mm10)
  a <- intersect(list_hg38,list_mm10) #intersect: 两个向量的交集，集合可以是数字、字符串等,返回对应的值
  
  if(length(a) != 0){
    hg38_tombo_0109_br <- readRDS('../file/hg38_tombo_0109_br.rds')
    hg38_tombo_0109[a]$brts_score <- hg38_tombo_0109_br[a]$br_score #对应上给分与br_score一样
  }
  
}

#-------------------------------------------------------------------------------
saveRDS(hg38_tombo_0109,'../file/hg38_tombo_0109_ts.rds')


