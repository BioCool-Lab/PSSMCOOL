#' Singular Value Decomposition (SVD)
#' @description Singular value decomposition is a general purpose matrix factorization approach
#'that has many useful applications in signal processing and statistics. In this function SVD is
#'applied to a matrix representation of a protein with the aim of reducing its dimensionality
#'Given an input matrix Mat with dimensions N*M SVD is used to calculate its factorization
#'of the form: \eqn{Mat=U\Sigma V,} where \eqn{\Sigma} is a diagonal matrix whose diagonal
#'entries are known as the singular values of Mat. The resulting descriptor is the ordered
#'set of singular values: \eqn{SVD\in\mathcal{R}^L,} where L=min(M,N).
#'and here \code{\link[base]{svd}} function is used for this purpose.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return feature vector of length 20
#' @references
#' L. Nanni, A. Lumini, and S. J. T. S. W. J. Brahnam, "An empirical study of different approaches for protein classification,"
#' vol. 2014, 2014.
#' @export
#' @examples
#' w<-SVD_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
SVD_PSSM<-function(pssm_name){
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
  SVD_vec<-svd(p)
  SVD_vec<-SVD_vec$d
  SVD_vec<-round(SVD_vec,digits = 3)
  return(SVD_vec)
}
