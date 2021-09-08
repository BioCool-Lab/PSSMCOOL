#' Linear predictive coding feature
#' @description This function uses Linear predictive coding algorithm for each column of PSSM Matrix
#'. so in this script \code{\link[phonTools]{lpc}} function is used which produces a 14-dimensional
#'vector for each column, since PSSM has 20 column eventually it will be obtained a 20*14=280 dimensional
#'feature vector for each PSSM Matrix by this function.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @importFrom phonTools lpc
#' @return feature vector of length 280
#' @references
#' L. Li et al., "PSSP-RFE: accurate prediction of protein structural class by recursive feature extraction
#' from PSI-BLAST profile, physical-chemical property and functional annotations," vol. 9, no. 3, 2014.
#' @export
#' @examples
#' w<-LPC_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
LPC_PSSM<-function(pssm_name){
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
  E<-x
  E<-1/(1+exp(-E))
  requireNamespace("phonTools",quietly = TRUE)
  LPC_vec<-rep(0,280)
  for(i in 1:20){
    LPC_vec[((i-1)*14+1):((i-1)*14+14)]<-phonTools::lpc(E[,i])
  }
  LPC_vec<-round(LPC_vec,digits = 4)
  return(LPC_vec)
}
