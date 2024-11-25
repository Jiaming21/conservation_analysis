#生成染色体信息
#hg19
  #定义GFgenreation_m6A函数
  #background是hg19InMm10_m6AConservationInfo_updated_3.rds的前2400行
  #newdata更新, 将background和newdata合并
  #genome(newdata) <- NA #genome 属性为空
  #newdata更新，取前26列
  #GF_positiveAll是利用GFgenreation_m6A()函数从newdata里取的特征信息
  #GF_positiveAll更新为dataframe格式
# colnames(mcols(background))
# [1] "PubmedID"                          "GSE"                               "Technique"                        
# [4] "Cell_line"                         "Treatment"                         "site_pos"                         
# [7] "reference_sequence"                "A_41bp"                            "NumberOfExperimentSupported"      
# [10] "ID"                                "liftoverToMm10"                    "site_pos_mm10"                    
# [13] "conserved_status"                  "PubmedID_mm10"                     "GSE_mm10"                         
# [16] "Technique_mm10"                    "Cell_line_mm10"                    "Treatment_mm10"                   
# [19] "NumberOfExperimentSupported_mm10"  "ID_mm10"                           "tissue_specific"                  
# [22] "br_score"                          "brts_score"                        "nodsHg19_score"                   
# [25] "nodsmm10_score"                    "sequence_score"|                   "PUlearningHg19_score"             
# [28] "PUlearningMm10_score"              "DP2"                               "MeRIPTissueSpecific_hg19"         
# [31] "MeRIPTissueSpecific_mm10"          "MeRIPts_score"                     "Gene"                             

# colnames(mcols(newdata))
# [1] "PubmedID"                         "GSE"                              "Technique"                       
# [4] "Cell_line"                        "Treatment"                        "site_pos"                        
# [7] "reference_sequence"               "A_41bp"                           "NumberOfExperimentSupported"     
# [10] "ID"                               "liftoverToMm10"                   "site_pos_mm10"                   
# [13] "conserved_status"                 "PubmedID_mm10"                    "GSE_mm10"                        
# [16] "Technique_mm10"                   "Cell_line_mm10"                   "Treatment_mm10"                  
# [19] "NumberOfExperimentSupported_mm10" "ID_mm10"                          "tissue_specific"                 
# [22] "br_score"                         "brts_score"                       "nodsHg19_score"                  
# [25] "nodsmm10_score"                   "sequence_score"     
  #newdata更新，行合并background(前)和newdata(后)，列以background为准，没有信息自动填充NA

  #GF_positiveAll里有NA的列都赋值为0
  #移除第30列“motif DRACH”, 都为1


