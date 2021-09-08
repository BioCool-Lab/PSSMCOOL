#' Discrete Cosin Transform Feature
#' @description To construct this feature vector, Two-Dimensional DCT algorithm has been used by applying
#'\code{\link[dtt]{dct}} function from dtt package which DCT stands for Discrete Cosin Transform.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @importFrom dtt dct
#' @return feature vector of length 400
#' @references
#' Wang, L., et al., Advancing the prediction accuracy of protein-protein interactions by utilizing
#' evolutionary information from position-specific scoring matrix and ensemble classifier. 2017. 418: p. 105-110.
#'
#' Y. Wang, Y. Ding, F. Guo, L. Wei, and J. J. P. o. Tang, "Improved detection of DNA-binding proteins
#' via compression technology on PSSM information," vol. 12, no. 9, 2017.
#'
#' @export
#' @examples
#' as<-Discrete_Cosine_Transform(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
Discrete_Cosine_Transform<- function(pssm_name){
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
  w<-c()
  colnames(p)<-NULL
  p<-1/(1+exp(-p))
  requireNamespace("dtt",quietly = TRUE)
  p<-dtt::dct(p)
  p<-round(p,digits = 4)
  for(i in 1:20){
    w<-c(w,p[i,])
  }
  return(w)
}

