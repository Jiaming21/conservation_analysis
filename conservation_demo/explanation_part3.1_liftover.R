#mapWithDatabase_explanation

#newdata: 新发现的m6A, 关于metadata: 

  #人的信息都一样: 
  #"PubmedID"
  #"GSE"
  #"Technique"
  #"Cell_line"
  #"Treatment"
  #以上均为'-'
  #"site_pos"                        
  #"reference_sequence"
  #"A_41bp"
  #"NumberOfExperimentSupported"     
  #"ID"
  #填充人的对应信息

#鼠的信息要分类: 

#数据中有可以转到鼠的

  #对于不可以转到鼠的行: 

    #"liftoverToMm10": 'NO'
    #"conserved_status": 'not conserved'
    #"site_pos_mm10"

    #"PubmedID_mm10"
    #"GSE_mm10"                        
    #"Technique_mm10"
    #"Cell_line_mm10"
    #"Treatment_mm10"                  
    #"NumberOfExperimentSupported_mm10"
    #"ID_mm10“
    #以上均为'-'

  #对于可以转到鼠的行: 

    #可以在保守鼠m6A数据库中找到

      #"liftoverToMm10": 'YES'
      #"conserved_status": 'conserved'
      #"site_pos_mm10" #填充从hg19_to_mm10来的对应信息

      #"PubmedID_mm10"
      #"GSE_mm10"                        
      #"Technique_mm10"
      #"Cell_line_mm10"
      #"Treatment_mm10"                  
      #"NumberOfExperimentSupported_mm10"
      #"ID_mm10“
      #以上填充从m6A_SB_mm10来的对应信息

    #不可以在保守鼠m6A数据库中找到

      #"liftoverToMm10": 'YES'
      #"conserved_status": 'not conserved'
      #"site_pos_mm10" #填充从hg19_to_mm10来的对应信息

      #"PubmedID_mm10"
      #"GSE_mm10"                        
      #"Technique_mm10"
      #"Cell_line_mm10"
      #"Treatment_mm10"                  
      #"NumberOfExperimentSupported_mm10"
      #"ID_mm10“
      #以上均为'-'

#数据中全部都不可以转到鼠

#"liftoverToMm10": 'NO'
#"conserved_status": 'not conserved'
#"site_pos_mm10"

#"PubmedID_mm10"
#"GSE_mm10"                        
#"Technique_mm10"
#"Cell_line_mm10"
#"Treatment_mm10"                  
#"NumberOfExperimentSupported_mm10"
#"ID_mm10“
#以上均为'-'


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

#site_pos_mm10: 存鼠的位点信息

#m6A_SB_mm10: 鼠保守m6A数据库
