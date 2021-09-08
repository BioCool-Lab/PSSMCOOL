#' smoothed PSSM feature
#' @description In this function at first a Matrix called smoothed-PSSM is constructed from PSSM Matrix
#'by applying "ws" parameter which called sliding window size and taken from user and usually is equals to 7. Then
#'using other window size parameter "w" which usually equals to 11 at each position smoothed feature vector is
#'constructed.
#' @param pssm_name name of PSSM Matrix file
#' @param ws window size for smoothing PSSM Matrix
#' @param w window size for extracting feature vector
#' @param v vector of desired positions to extract their features
#' @import utils
#' @details In the construction of a smoothed PSSM, each row vector of a residue \eqn{\alpha_i} is represented and
#'smoothed by the summation of ws surrounding row vectors \eqn{(V_{smoothed_i}=V_{i-(ws-1)/2}+...+V_i+...+
#'V_{i+(ws+1)/2})} For the N-terminal and C-terminal of a protein, (w-1)/2 ZERO vectors, are appended to the
#'head or tail of a smoothed PSSM profile. Using the smoothed PSSM encoding scheme the feature vector of a residue
#'\eqn{\alpha_i} is represented by \eqn{(V_{smoothed_i-(ws-1)/2},...,V_{smoothed_i},...,V_{smoothed_i+(ws+1)/2})}
#'The feature values in each vector are normalized to a range between -1 and 1.
#' @return a matrix of feature vectors
#' @references
#' Cheng, C.W., et al. (2008) Predicting RNA-binding sites of proteins using support vector machines and
#' evolutionary information, BMC Bioinformatics, 9 Suppl 12, S6.
#' @seealso \code{\link{kiderafactor}}
#' @export
#' @examples
#' w<-smoothed_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7,11,c(2,3,8,9))
smoothed_PSSM<-function(pssm_name,ws,w,v=NULL){ #ws in range 3,5,...,11 and w is in range 3,5,...,41
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
  m<-x
  m<-1/(1+exp(-m))
  L<-dim(m)[1]
  smoothed_PSSM<-matrix(0,L,20)
  h<-matrix(0,(ws-1)/2,20)
  k<-matrix(0,(ws-1)/2,20)
  E<-c()
  m<-rbind(h,m,k)
  for(i in 1:L){
    for(j in 0:(ws-1)){
      E<-rbind(E,m[i+j,])
    }
    smoothed_PSSM[i,]<-colSums(E)
    E<-c()
  }
  mh<-matrix(0,(w-1)/2,20)
  mk<-matrix(0,(w-1)/2,20)
  M<-rbind(mh,smoothed_PSSM,mk)
  d<-(w-1)/2 +1
  a<-1
  x<-matrix(0,L,20*w)
  w1=w2<-c()
  for(i in (d-(w-1)/2):(d+(w-1)/2)){
    w2<-c(w2,M[i,])
  }
  x[a,]<-w2
  a<-a+1
  d<-i+1
  while(d<=(L+(w-1))){
    w2<-c(w2[21:length(w2)],M[d,])
    x[a,]<-w2
    a<-a+1
    d<-d+1

  }
  if(length(v)!=0){
    y<-x[v,]
    y<-as.data.frame(y)
    rownames(y)<-v
    colnames(y)<-1:(w*20)
  }
  else{
    y<-x
    y<-as.data.frame(y)
    rownames(y)<-1:L
    colnames(y)<-1:(w*20)
  }
  return(round(y,digits = 4))
}
