#上接part3.4_numberOfDataset_explanation: 
#newdata有了"nodsHg19_score"和"nodsmm10_score"
#hg19_to_mm10有了"br_score"；完美匹配，偏移匹配的行在此步去除；只剩下可以转换到鼠但无法匹配的行
#---------------------------------------------------------------------------------------------#
#newdata: 新发现的m6A, 关于metadata: 
#人的信息: 
#"PubmedID"
#"GSE"
#"Technique"
#"Cell_line"
#"Treatment"
#以上填充'-'
#"site_pos"                        
#"reference_sequence"
#"A_41bp"
#"NumberOfExperimentSupported"     1
#"ID"
#以上填充人的对应信息
                                           #偏移匹配不上     #偏移匹配上了     #完美匹配
#"liftoverToMm10":       'NO'              'YES'             'YES'             'YES'
#"conserved_status":     'not conserved'   'not conserved'   'not conserved'   'conserved' 
#"site_pos_mm10":        '-'               '鼠的位点'        '鼠的位点'        '鼠的位点'

#"PubmedID_mm10"         '-'               '-'               '鼠的信息'        '鼠的信息'
#"GSE_mm10"              '-'               '-'               '鼠的信息'        '鼠的信息'         
#"Technique_mm10"        '-'               '-'               '鼠的信息'        '鼠的信息'
#"Cell_line_mm10"        '-'               '-'               '鼠的信息'        '鼠的信息'
#"Treatment_mm10"        '-'               '-'               '鼠的信息'        '鼠的信息'
#"NumberOfExperimentSupported_mm10" 1      1                 '鼠的信息'        '鼠的信息'
#"ID_mm10“               '-'               “挑选”过的人ID   “挑选”过的人ID   “挑选”过的人ID
#"br_score"                   0                 0            [0.8,0.6,0.4,0.2]     1
#"brts_score"                 0                 0       部分[0.8,0.6,0.4,0.2] 部分[0.8,0.6,0.4,0.2]
#以上如有信息，来自于m6A_SB_mm10保守库
#"nodsHg19_score"       （文献数量）全部为1，所以全部给0.2 
#"nodsmm10_score"       （文献数量）默认为0，1给0.2，2给0.4，3-4给0.6，5-6给0.8，大于6给1

#---------------------------------------------------------------------------------------------#
#上接part3.2_base-resolution_explanation
#hg19_to_mm10: 经转换得到的鼠的位点, metadata(都是人)有: 
#"PubmedID"(人)
#"GSE"(人)
#"Technique"(人)
#"Cell_line"(人)
#"Treatment"(人)
#以上均为'-'
#"site_pos"(人鼠一样，不符合的人map后没了，一直做的是删减工作)                        
#"reference_sequence"(人)       
#"A_41bp"(人)  
#"NumberOfExperimentSupported"(人)      
#"ID"(人鼠一样，不符合的人map后没了，一直做的是删减工作)
#site_pos_mm10(鼠)
#"br_score": 0，【0.8，0.6，0.4，0.2】偏移符合行
#完美匹配，偏移匹配的行在此步去除；只剩下可以转换到鼠但无法匹配的行
#-------------------------------------------------------------------#
#hg19_to_mm10: 鼠的位点(全部碱基都是A), 关于metadata: 
#人的信息: 
#"PubmedID"
#"GSE"
#"Technique"
#"Cell_line"
#"Treatment"
#以上填充'-'
#"site_pos"                        
#"reference_sequence"
#"A_41bp"
#"NumberOfExperimentSupported"     1
#"ID"
#以上填充人的对应信息
                                           #偏移匹配不上     #偏移匹配上了     #完美匹配
#"liftoverToMm10":                            'YES'             'YES'             'YES'
#"conserved_status":                     'not conserved'      'conserved'      'conserved' 
#"site_pos_mm10":                           '鼠的位点'        '鼠的位点'        '鼠的位点'
 
#"PubmedID_mm10"                            '-'               '鼠的信息'        '鼠的信息'
#"GSE_mm10"                                 '-'               '鼠的信息'        '鼠的信息'         
#"Technique_mm10"                           '-'               '鼠的信息'        '鼠的信息'
#"Cell_line_mm10"                           '-'               '鼠的信息'        '鼠的信息'
#"Treatment_mm10"                           '-'               '鼠的信息'        '鼠的信息'
#"NumberOfExperimentSupported_mm10"         1                 '鼠的信息'        '鼠的信息'
#"ID_mm10“                               “挑选”过的人ID   “挑选”过的人ID     “挑选”过的人ID
#"br_score"                                 0              [0.8,0.6,0.4,0.2]        1
#"brts_score"                                  0       部分[0.8,0.6,0.4,0.2] 部分[0.8,0.6,0.4,0.2]
#以上如有信息，来自于m6A_SB_mm10保守库
#"nodsHg19_score"       （文献数量）全部为1，所以全部给0.2 
#"nodsmm10_score"       （文献数量）默认为0，1给0.2，2给0.4，3-4给0.6，5-6给0.8，大于6给1
#"A_11bp"      #人的对应位置的11bp的序列
#"A_mm10"      #鼠碱基
  #删去"A_mm10"不是A的行
#"A_11bp_mm10" #鼠的对应位置的11bp的序列
#"status"                               'not conserved'      'not conserved'   'conserved'
#"score"：此时hg19_to_mm10全部碱基都是A，全部都有分

#-------------------------------------------------------------------#
#A_11bp("DNAStringSet"): newdata对应的人11bp序列
#hg19_to_mm10: 新的鼠(由newdata chain转换而来)

#有可以转换到鼠的
#hg19_to_mm10$ID newdata$ID 都是人的信息，hg19_to_mm10$ID只保留了转过去的
#match: 前者在后者的行(后者完全包括前者)
#hg19_to_mm10$A_11bp: match到的人11bp序列
#hg19_to_mm10$A_mm10: 根据hg19_to_mm10坐标转换到的鼠碱基(为了确定是否是A)
#indx: hg19_to_mm10中对应位置是A的行
  #有转过去是A的
  #hg19_to_mm10只剩转过去是A的行
  #hg19_to_mm10$A_11bp_mm10，鼠的对应位置的11bp的序列
  #hg19_to_mm10$status默认为'-'，与鼠保守库match取出保守行用'conserved'
  
  #hg19_to_mm10$score默认为0
  #SM_motif  #SM_motif 区域的打分矩阵
  #  A G C T
  #A 4 3 0 0
  #G 3 4 0 0
  #C 0 0 4 3
  #T 0 0 3 4

  #SM        #其他区域的打分矩阵
  #  A G C T
  #A 2 1 0 0
  #G 1 2 0 0
  #C 0 0 2 1
  #T 0 0 1 2

  #for
  #seq_hg19 <- hg19_to_mm10$A_11bp[[i]]       #[[i]]取出11bp的人序列
  #seq_mm10 <- hg19_to_mm10$A_11bp_mm10[[i]]  #[[i]]取出11bp的鼠序列
    #for
    #每一个碱基打分，总分为score_1 
    #score_normalized为score_1/32 #归一化，除以最大分数
    #hg19_to_mm10$score[i] 对应行给上归一化的分数
    #显示打分进度

  #newdata$sequence_score默认为0
  #将newdata$ID match hg19_to_mm10$ID的行给上hg19_to_mm10的score

  #转过去全都不是A
  #newdata$sequence_score为0
#都不可以转换到鼠
#newdata$sequence_score为0

#hg19_to_mm10
#从转换到的 到 碱基是A的

#现在的newdata：
#newdata: 新发现的m6A, 关于metadata: 
#人的信息: 
#"PubmedID"
#"GSE"
#"Technique"
#"Cell_line"
#"Treatment"
#以上填充'-'
#"site_pos"                        
#"reference_sequence"
#"A_41bp"
#"NumberOfExperimentSupported"     1
#"ID"
#以上填充人的对应信息
#偏移匹配不上     #偏移匹配上了     #完美匹配
#"liftoverToMm10":       'NO'              'YES'             'YES'             'YES'
#"conserved_status":     'not conserved'   'not conserved'   'not conserved'   'conserved' 
#"site_pos_mm10":        '-'               '鼠的位点'        '鼠的位点'        '鼠的位点'

#"PubmedID_mm10"         '-'               '-'               '鼠的信息'        '鼠的信息'
#"GSE_mm10"              '-'               '-'               '鼠的信息'        '鼠的信息'         
#"Technique_mm10"        '-'               '-'               '鼠的信息'        '鼠的信息'
#"Cell_line_mm10"        '-'               '-'               '鼠的信息'        '鼠的信息'
#"Treatment_mm10"        '-'               '-'               '鼠的信息'        '鼠的信息'
#"NumberOfExperimentSupported_mm10" 1      1                 '鼠的信息'        '鼠的信息'
#"ID_mm10“               '-'               “挑选”过的人ID   “挑选”过的人ID   “挑选”过的人ID
#"br_score"                   0                 0            [0.8,0.6,0.4,0.2]     1
#"brts_score"                 0                 0       部分[0.8,0.6,0.4,0.2] 部分[0.8,0.6,0.4,0.2]
#以上如有信息，来自于m6A_SB_mm10保守库
#"nodsHg19_score"       （文献数量）全部为1，所以全部给0.2 
#"nodsmm10_score"       （文献数量）默认为0，1给0.2，2给0.4，3-4给0.6，5-6给0.8，大于6给1
#"sequence_score"             0                 0/0-1           0/0-1            0/0-1
















