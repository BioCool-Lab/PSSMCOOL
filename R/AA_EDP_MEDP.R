#' EDP_EEDP_MEDP feature vector
#' @description This is a feature vector of length 420 which is used for prediction of protein
#'structural class for low-similarity sequences.at first ED-PSSM Matrix with 20*20 dimensions
#'is constructed from PSSM Matrix then by using this Matrix, EDP and EEDP vectors are
#'obtained eventually MEDP feature vector is obtained by fusing these vectors.
#' @param pssm_name is name of PSSM Matrix  file
#' @import utils
#' @return a feature vectors of length 420
#' @references
#' Zhang, L., Zhao, X. and Kong, L. (2014) Predict protein structural class for low-similarity sequences by evolutionary difference
#' information into the general form of Chou's pseudo amino acid composition, Journal of Theoretical Biology, 355, 105-110.
#'
#' @export
#'
#' @examples
#' as<-EDP_MEDP(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GS61.txt.pssm"))
EDP_MEDP <- function(pssm_name){
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
  p<-x
  p<-1/(1+exp(-p))
  L<-dim(p)[1]
  s<-0
  e<-matrix(0,20,20)
  for(k in 1:20){
    for(t in 1:20){
      for(i in 2:(L-1)){
        edf<-(p[i-1,k]-p[i+1,t])/2
        edf<-edf^2
        s<-s+edf

      }
      e[k,t]<-s/(L-2)
      s<-0
    }
  }
  EDP<-apply(e,2,mean)
  names(EDP)<-NULL
  EDP<-round(EDP,digits = 4)
  v<-c()
  for(i in 1:20){
    v<-c(v,e[i,])
  }
  EEDP<-round(v,digits = 4)
  MEDP<-c(EDP,EEDP)
  return(MEDP)
}

