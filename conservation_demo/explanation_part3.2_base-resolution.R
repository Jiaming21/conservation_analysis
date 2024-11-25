#上接liftover-explanation

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
#"NumberOfExperimentSupported"     
#"ID"
#以上填充人的对应信息

#"liftoverToMm10":       'NO'              'YES'             'YES'
#"conserved_status":     'not conserved'   'not conserved'   'conserved'
#"site_pos_mm10":        '-'               '鼠的位点'        '鼠的位点'

#"PubmedID_mm10"         '-'               '-'               '鼠的信息'
#"GSE_mm10"              '-'               '-'               '鼠的信息'              
#"Technique_mm10"        '-'               '-'               '鼠的信息'
#"Cell_line_mm10"        '-'               '-'               '鼠的信息'
#"Treatment_mm10"        '-'               '-'               '鼠的信息'
#"NumberOfExperimentSupported_mm10" 1      1                 1
#"ID_mm10“               '-'               “挑选”过的人num   “挑选”过的人num  
#以上如有信息，来自于m6A_SB_mm10保守库


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


#m6A_SB_mm10: 鼠保守库
#"PubmedID"
#"GSE"
#"Technique"
#"Cell_line"
#"Treatment"
#"site_pos"                        
#"reference_sequence"     
#"A_41bp"
#"NumberOfExperimentSupported"    
#"ID"
#以上都填充鼠的信息

#经过base-resolution变成:

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
#"NumberOfExperimentSupported"     
#"ID"
#以上填充人的对应信息

#"liftoverToMm10":       'NO'              'YES'             'YES'
#"conserved_status":     'not conserved'   'not conserved'   'conserved'
#"site_pos_mm10":        '-'               '鼠的位点'        '鼠的位点'

#"PubmedID_mm10"         '-'               '-'               '鼠的信息'
#"GSE_mm10"              '-'               '-'               '鼠的信息'              
#"Technique_mm10"        '-'               '-'               '鼠的信息'
#"Cell_line_mm10"        '-'               '-'               '鼠的信息'
#"Treatment_mm10"        '-'               '-'               '鼠的信息'
#"NumberOfExperimentSupported_mm10" 1      1                 '鼠的信息'
#"ID_mm10“               '-'               “挑选”过的人num   “挑选”过的人num
#"br_score"              0                 0                 1
#以上如有信息，来自于m6A_SB_mm10保守库


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

#保守行删去(hg19_to_mm10的site_pos_mm10(鼠)与保守库site_pos匹配)

#newdata加上了“新保守行”的信息
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
#"NumberOfExperimentSupported"     
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
#以上如有信息，来自于m6A_SB_mm10保守库






