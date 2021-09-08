#' DPC-PSSM,AAC-PSSM and AADP-PSSM feature vectors
#' @description This feature is combination of amino asid composition and dipeptide composition feature
#'vectors. DPC feature stands for dipeptide composition, which multiplies the values in two consecutive
#'rows and two different columns, calculating this for each of the different columns and obtaining the
#'sum of these and for each. And for both columns, the product divides their sum by L-1, and because the
#'result depends on two different columns, length of this feature vector would be 400. AAC-PSSM is actually
#'mean of PSSM Matrix columns which its length is 20. eventually AADP-PSSM is combination of these vectors
#'and with length 420.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#'
#' @return feature vector of length 420
#' @references
#' Liu, T., Zheng, X. and Wang, J. (2010) Prediction of protein structural class for
#' low-similarity sequences using support vector machine and PSI-BLAST profile, Biochimie, 92, 1330-1334.
#' @export
#' @examples
#' ss<-DPC_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
DPC_PSSM <- function(pssm_name){
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
  AAC<-apply(m2,2,mean)
  names(AAC)<-NULL
  AAC<-round(AAC,digits = 4)
  y<-matrix(0,nrow = 20,ncol = 20)
  for (i in 1:20) {
    for (j in 1:20) {
      for (k in 1:(L-1)) {
        y[i,j]<-y[i,j]+m2[k,i]*m2[k+1,j]
      }
    }
  }
  y<-(1/(L-1))*y
  vec<-c()
  for(i in 1:20){
    vec<-c(vec,y[i,])
  }
  vec<-round(vec,digits = 4)
  AADP<-c(AAC,vec)
  return(AADP)
}
