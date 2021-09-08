#' k_seperated_bigrame feature vector
#' @description This feature is almost identical to the \code{\link{DPC_PSSM}} feature, and in fact the DPC feature is part of
#'this feature (for k=1) and for two different columns, considers rows that differ by the size of the unit k.
#' @param pssm_name is name of PSSM Matrix file
#' @param k a parameter that specifies separated length between amino acids
#' @import utils
#' @return a feature vector of length 400
#' @references
#' Saini, H., et al.(2016) Protein Fold Recognition Using Genetic Algorithm Optimized Voting Scheme and Profile
#' Bigram.
#' @export
#' @examples
#'  w<-k_seperated_bigrame(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),5)
k_seperated_bigrame<-function(pssm_name,k){
  x<-read.delim(pssm_name,skip = 2,sep = "",header = FALSE)
  x<-x[-1,-c(1,23:44)]
  d<-which(x=="Lambda")
  if(length(d)!=0){
    x<-x[-c(d:dim(x)[1]),]
  }
  x<-x[,-1]
  colnames(x)<-NULL
  rownames(x)<-NULL
  x<-as.matrix(x)
  mode(x)<-"integer"
  m2<-x
  m2<-1/(1+exp(-m2))
  L<-dim(m2)[1]
  j<-1
  km<-matrix(0,20,20)
  s<-0
    for (m in 1:20) {
      for (n in 1:20) {
        for (i in 1:(L-k)) {
          s <- s+m2[i,m]*m2[i+k,n]
        }
        km[m,n]<-s
        s<-0
      }
    }
  v<-c()
  for(i in 1:20){
    v<-c(v,km[i,])
  }
  v<-round(v,digits = 4)
  return(v)
}




