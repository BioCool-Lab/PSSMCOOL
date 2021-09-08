#' 3-Mer and 2-Mer
#' @description This function produces all possible k-mers from 20 amino acids for use in other functions.
#' @param k is length of k-mer which user imports
#' @import utils
#' @import gtools
#' @return a matrix which its first row includes all k-mers
#' @export
#'
#' @examples
#' ax<-three_mer(3)
three_mer<-function(k){
  v<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  requireNamespace("gtools",quietly = TRUE)
  if(k==3){
      aq<-gtools::permutations(20,3,v,repeats.allowed = T)
      ww<-t(aq)
      aa<-apply(ww,2,function(y) paste0(y,collapse = ""))
      az<-rep(0,8000)
      as<-rbind(aa,az)
      rownames(as)=NULL
  } else if(k==2){
      aq<-gtools::permutations(20,2,v,repeats.allowed = T)
      ww<-t(aq)
      aa<-apply(ww,2,function(y) paste0(y,collapse = ""))
      az<-rep(0,400)
      as<-rbind(aa,az)
      rownames(as)=NULL
  }
  else{
    stop("ERROR: enterd Number must be 2 or 3")
  }
  return(as)
}
############################################################################################################
#' 3-mer and 2-mer in dataframe
#' @description This function produces all possible 2-mers or 3-mers by counting paths of length 2 or 3 in
#'a dataframe which is thought as a graph
#' @param s a dataframe with 2 columns
#' @param h is length of k-mer
#' @import utils
#' @return all k-mers by counting paths of length h in dataframe which is considered as a graph
#' @export
#'
#' @examples
#' s1<-LETTERS[1:4]
#' s2<-LETTERS[3:6]
#' s<-data.frame(s1,s2)
#' dc<-k_mers(s,3)
k_mers<-function(s,h){
  aq<-c()
  if(h==3){
    for(i in 1:(dim(s)[1]-h+1)){
      aq<-c(aq,paste0(s[i,1],s[i+1,1],s[i+2,1]),paste0(s[i,1],s[i+1,1],s[i+2,2]),paste0(s[i,1],s[i+1,2],s[i+2,2]),paste0(s[i,1],s[i+1,2],s[i+2,1]))
      aq<-c(aq,paste0(s[i,2],s[i+1,2],s[i+2,2]),paste0(s[i,2],s[i+1,2],s[i+2,1]),paste0(s[i,2],s[i+1,1],s[i+2,2]),paste0(s[i,2],s[i+1,1],s[i+2,1]))
    }
  }
  else if(h==2){
    for(i in 1:(dim(s)[1]-h+1)){
      aq<-c(aq,paste0(s[i,1],s[i+1,1]),paste0(s[i,1],s[i+1,2]),paste0(s[i,2],s[i+1,2]),paste0(s[i,2],s[i+1,1]))
    }
  }
  else{
    stop("ERROR: Enterd Number must be 2 or 3")
  }
  return(aq)
}
######################################################################################
#' SCSH Feature vector
#' @description This function gets a PSSM Matrix as input and extracts corresponding protein and consunsus
#'sequence from it. By placing these two vectors next to each other, a dataframe is created.
#'In each row of this dataframe, each component is connected to the next row components by an arrow. so a
#'directed graph is produced. next by use of previous functions a feature vector of length 400 or 8000 is created.
#'
#' @param pssm_name is name of PSSM Matrix file
#' @param k a parameter indicates length of k-mer
#' @import utils
#' @return feature vector of length 400 or 8000
#' @references
#' Zahiri, J., et al., LocFuse: human protein–protein interaction prediction via classifier fusion using protein
#' localization information. Genomics, 2014. 104(6): p. 496-503.
#' @export
#'
#' @examples
#' zz<- scsh2(system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL"),2)
scsh2<-function(pssm_name,k){
  x<-read.delim(pssm_name,skip = 2,sep = "",header = FALSE)
  x<-x[-1,-c(1,23:44)]
  d<-which(x=="Lambda")
  if(length(d)!=0){
    x<-x[-c(d:dim(x)[1]),]
  }
  colnames(x)<-NULL
  rownames(x)<-NULL
  x<-as.matrix(x)
  k2<-x[,1]
  k2<-as.character(k2)
  L<-length(k2)
  w<-consunsus_sequence(pssm_name)
  w2<-unlist(strsplit(w,""))
  s<-data.frame(k2,w2)
  kmer<-k_mers(s,k)
  mk<-three_mer(k)
  for(i in 1:length(kmer)){
    kl<-which(mk[1,]==kmer[i])
    mk[2,kl]<-1
  }
  vv<-as.integer(mk[2,])
  return(vv)

}

