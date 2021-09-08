#' AATP_TPC feature vector
#' @description For getting this feature which was used to protein structural class prediction,
#'at first mean of every column in PSSM Matrix is computed to achieve a 20-dimensional vector
#'called AAC.then by fusing it with other vector of length 400 called TPC, which is similar to \code{\link{DPC_PSSM}}
#'AATP feature vector of length 420 is obtained.
#' @param pssm_name is name of PSSM Matrix file
#' @import utils
#' @return a feature vector of length 420
#' @references
#' Zhang, S., Ye, F. and Yuan, X. (2012) Using principal component analysis and support vector machine to predict protein
#' structural class for low-similarity sequences via PSSM, Journal of Biomolecular Structure & Dynamics, 29, 634-642.
#' @seealso \code{\link{DPC_PSSM}}
#' @export
#'
#' @examples
#' as<-AATP_TPCC(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GQS7.txt.pssm"))
AATP_TPCC <- function(pssm_name){
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
  X<-apply(p,2,mean)
  names(X)<-NULL
  y<-matrix(0,20,20)
  AATP<-c()
  for(i in 1:20){
    kjey<-rep(0,20)
    for(j in 1:20){
      for(k in 1:(L-1)){
        kjey[j]<-kjey[j]+p[k,i]*p[k+1,j]
      }
    }
    for(j in 1:20){
      y[i,j]<-kjey[j]/sum(kjey)
      AATP<-c(AATP,y[i,j])
    }
  }
  AATPTPc<-c(X,AATP)
  AATPTPc<-round(AATPTPc,digits = 4)
  return(AATPTPc)
}
