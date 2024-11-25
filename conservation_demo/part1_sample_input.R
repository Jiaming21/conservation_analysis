#前情回顾: sample_input是输入样本
#Convert to Granges，将输入转成Granges格式
PIANO_INPUT <- GRanges(seqnames = paste0('chr',as.character(sample_input$V1)),
                                    IRanges(start = as.numeric(sample_input$V2),
                                            end = as.numeric(sample_input$V3)),
                                    strand = sample_input$V4)

PIANO_INPUT <- m6ALogisticModel::sample_sequence("A",PIANO_INPUT,Hsapiens,Fixed = F) #处理输入文件中的区间，查找Granges范围内的等于碱基A的位置并进行输出
#可以查看length(PIANO_INPUT)=8072; unique(PIANO_INPUT@ranges@width)=1(单个碱基)

PIANO_INPUT <- m6ALogisticModel::sample_sequence("DRACH",PIANO_INPUT + 2,Hsapiens,Fixed = F)  #m6A甲基化修饰主要受METT3/METTL14酶调控，修饰特征序列为DRACH（D:A\G\U，R=A\G，H=A\C\U）
#但是注意，并不是所有符合该序列规则的序列都能被甲基化，仅有小部分含有DRACH motif序列的片段能发生m6A修饰
#length(PIANO_INPUT)=546

PIANO_INPUT <- PIANO_INPUT - 2 #PIANO_INPUT@ranges
print(paste0('length of PIANO_INPUT: ',length(PIANO_INPUT)))
saveRDS(PIANO_INPUT,'./PIANO_INPUT.rds')
rm(sample_input)