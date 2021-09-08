#' discrete wavelet transform feature vector
#' @description In construction of this feature vector, the \code{\link[waveslim]{dwt.nondyadic}} function is used from
#'"waveslim", package
#'to calculate the discrete wavelet transform for each column of the PSSM matrix, which considers it as a
#'discrete signal. At last, 4 levels DWT is used to analysis of these discrete signals of PSSM (each column) and
#'extracted the PSSM-DWT feature from PSSM of protein.
#' @param pssm_name name of pssm Matrix file
#' @import utils
#' @return feature vector of length 80
#' @references
#' Y. Wang, Y. Ding, F. Guo, L. Wei, and J. J. P. o. Tang, "Improved detection of DNA-binding proteins
#' via compression technology on PSSM information," vol. 12, no. 9, 2017.
#'
#' Y. Wang, Y. Ding, J. Tang, Y. Dai, F. J. I. A. t. o. c. b. Guo, and bioinformatics, "CrystalM: a
#' multi-view fusion approach for protein crystallization prediction," 2019.
#'
#' @export
#' @examples
#' as<-dwt_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
dwt_PSSM<-function(pssm_name){
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
  L<-dim(E)[1]
  requireNamespace("waveslim",quietly = TRUE)
  feature_vect<-c()
  for(j in 1:20){
    d<-c()
    a<-waveslim::dwt.nondyadic(E[,j])
    d<-c(d,a$d1,a$d2,a$d3,a$d4,a$s4)
    d<-round(d,digits = 4)
    s<-c(min(d),max(d),mean(d),sd(d))
    s<-round(s,digits = 4)
    feature_vect<-c(feature_vect,s)
  }
  return(feature_vect)
}

