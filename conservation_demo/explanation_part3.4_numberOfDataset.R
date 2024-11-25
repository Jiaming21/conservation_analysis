#上接tissue-specific_explanation: newdata有了brts_score
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
#"brts_score"                 0                 0       部分[0.8,0.6,0.4,0.2] 部分[0.8,0.6,0.4,0.2]
#以上如有信息，来自于m6A_SB_mm10保守库

#---------------------------------------------------------------------------------------------#
#在此步中，newdata中既有人的信息，也有鼠的信息；依照人的文献数量和鼠的文献数量给予"nodsHg19_score"和"nodsmm10_score"
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




