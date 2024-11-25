#前情回顾: 
#PIANO_INPUT: 符合“DRACH”的单个碱基位置, 由sample_input得来

hg19InMm10_m6AConservationInfo<-readRDS('./hg19InMm10_m6AConservationInfo_updated_3.rds') #已知的保守位点在人当中的数据库

hg19InMm10_m6AConservationInfo <- hg19InMm10_m6AConservationInfo[,-c(30:99)] 

database_overlap <- hg19InMm10_m6AConservationInfo[queryHits(findOverlaps(hg19InMm10_m6AConservationInfo,PIANO_INPUT))] 
# 输出值均为行号，根据行号取子集                           
# findOverlaps: 找交集
# queryHits: query中的交集数据的index; subjectHits: subject中的交集数据的index

newdata <- PIANO_INPUT[-queryHits(findOverlaps(PIANO_INPUT,database_overlap))] #newdata: 新发现的m6A
#取出PIANO_INPUT中不存在hg19InMm10_m6AConservationInfo的位置信息，invert = T, 反向。
#代码作用同: newdata <- subsetByOverlaps(PIANO_INPUT,database_overlap,invert = T)

print(paste0('length of newdata: ',length(newdata)))

print(paste0('length of database_overlap: ',length(database_overlap)))

rm(hg19InMm10_m6AConservationInfo)

rm(database_overlap)

rm(PIANO_INPUT)

