library(BSgenome)

NC1 <- function(data){
  for(i in 1:length(data) ){
    if(i==1){
      DCfirst <- unlist(as.vector(strsplit(data[1],"",fixed = TRUE)))
      DCsecond <- matrix(NA,nrow = length(data),ncol = length(DCfirst))
      DCsecond[1,] <-  DCfirst 
    }else{
      DCsecond[i,] <- unlist(as.vector(strsplit(data[i],"",fixed = TRUE)))
    }
  }
  return(DCsecond)
}

library(matrixStats)


CONPOSITION <- function(Data,NI=3,NTYPE="RNA",Freq=2){
  # MA <- NC1(DATA) 
  if(NTYPE=="RNA"){
    U="U"
  }else{
    U="T"
  }
  for(i in 1:NI){
    if(i==1){
      NP <- c("A","G","C",U)
    }else{
      NP1 <- NULL
      for(j in c("A","G","C",U)){
        for(k in 1:length(NP)){
          NP1 <- c(NP1,paste0(NP[k],j))
        }
      }
      NP <- NP1
    }
  }
  MA2 <- matrix(NA,ncol = length(NP),nrow = length(Data))
  colnames(MA2) <- NP
  if(NTYPE=="RNA"){
    for(i in 1:length(NP)){
      MA2[,i]<- vcountPattern(NP[i],RNAStringSet(Data))
    }
  }else{
    for(i in 1:length(NP)){
      MA2[,i]<- vcountPattern(NP[i],DNAStringSet(Data))
    }
  }
  M3 <- MA2
  if(Freq==1){
    for(i in 1:nrow(M3)){
      M3[i,] <-  MA2[i,]/sum(MA2[i,])
    }
  }
  return(M3)
}


sequenceFeatures <- function(Data,NTYPE="RNA") {
  sequences_M <- NC1(Data)
  if(NTYPE=="RNA"){
    U="U"
  }else{
    U="T"
  }
  N = ncol(sequences_M)
  
  cumFreq_A <- rowCumsums(matrix(as.numeric( sequences_M == "A" ), ncol = N, byrow = F))
  cumFreq_T <- rowCumsums(matrix(as.numeric( sequences_M ==  U ), ncol = N, byrow = F))
  cumFreq_C <- rowCumsums(matrix(as.numeric( sequences_M == "C" ), ncol = N, byrow = F))
  cumFreq_G <- rowCumsums(matrix(as.numeric( sequences_M == "G" ), ncol = N, byrow = F))
  
  cumFreq_combined <- matrix(0,ncol = N, nrow = length(Data))
  cumFreq_combined[sequences_M == "A"] <- cumFreq_A[sequences_M == "A"]
  cumFreq_combined[sequences_M == U] <- cumFreq_T[sequences_M == U]
  cumFreq_combined[sequences_M == "C"] <- cumFreq_C[sequences_M == "C"]
  cumFreq_combined[sequences_M == "G"] <- cumFreq_G[sequences_M == "G"]
  
  cumFreq_combined <- t(t(cumFreq_combined) / seq_len(N))
  colnames(cumFreq_combined ) <-  paste0("cumFreq_",seq_len(N))
  return(cumFreq_combined)
}


# TEST <- c("ACGUCUCUAUCGUACGUACGUAGUG","CUUCUCGUAACGUAGCAUACGGAUG")
# sequenceFeatures(TEST,"RNA")
# CONPOSITION(TEST)
