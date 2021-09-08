#' Averag Block feature vector
#' @description In this feature at first PSSM Matrix is divided to 20 Blocks. Then for each
#'Block mean of columns is computed to get 20-dimensional vector, eventually by appending these vectors to each other final feature
#'vector is obtained which would be of length 400. this feature vector is similar to \code{\link{PSSMBLOCK}} for N=20.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return feature vector of length 400
#' @references
#' L. Nanni, A. Lumini, and S. J. T. S. W. J. Brahnam, "An empirical study of different approaches for protein classification,"
#' vol. 2014, 2014.
#' @export
#' @examples
#' v<-Averag_Block(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
Averag_Block<-function(pssm_name){
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
  M<-floor(L/20)
  w<-c()
  for(i in 1:19){
    dd<-p[(1+(i-1)*M):(M+(i-1)*M),]
    v<-apply(dd,2,mean)
    names(v)<-NULL
    w<-c(w,v)
  }
  dd<-p[(1+19*M):L,]
  v<-apply(dd,2,mean)
  names(v)<-NULL
  w<-c(w,v)
  w<-round(w,digits = 4)
  return(w)
}
##################################################################################################
#' PSSM BLOCK feature vector
#' @description In this feature at first PSSM Matrix is divided to Blocks based on Number N which user imports. Then for each
#'Block mean of columns is computed to get 20-dimensional vector, eventually by appending these vectors to each other final feature
#'vector is obtained.
#' @param pssm_name neme of PSSM Matrix file
#' @param N number of blocks
#' @import utils
#' @return feature vector that it's length depends on parameter N
#' @references
#' J.-Y. An, L. Zhang, Y. Zhou, Y.-J. Zhao, and D.-F. J. J. o. c. Wang, "Computational methods using weighed-extreme learning
#' machine to predict protein self-interactions with protein evolutionary information," vol. 9, no. 1, p. 47, 2017.
#' @export
#' @examples
#' as<-PSSMBLOCK(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),5)
PSSMBLOCK <- function(pssm_name,N){
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
  w<-c();i<-1
  for(t in 1:N){
    M<-floor(L/t)
    while (i <=t){
      dd<-p[(1+(i-1)*M):(M+(i-1)*M),]
      v<-apply(dd,2,mean)
      names(v)<-NULL
      w<-c(w,v)
      i<-(i+1)
    }
    i<-1
  }
  w<-round(w,digits = 4)
  return(w)
}

