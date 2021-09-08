#' RPSSM feature
#' @description To obtain this feature, first the columns of the PSSM matrix are merged to obtain an L*10 matrix. Then,
#'with a relationship similar to the auto covariance transformation feature, this feature with a length of 110 is
#'obtained from this matrix.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return feature vector of length 110
#' @references
#' Ding, S., et al. (2014) A protein structural classes prediction method based on predicted secondary structure and
#' PSI-BLAST profile, Biochimie, 97, 60-65.
#' @export
#' @examples
#' w<-rpssm(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
rpssm<-function(pssm_name){
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
  L<-dim(m2)[1]
  m2<-1/(1+exp(-m2))
  p<-matrix(0,L,10)
  p[,1]<-(m2[,14]+m2[,19]+m2[,18])/3
  p[,2]<-(m2[,13]+m2[,11])/2
  p[,3]<-(m2[,10]+m2[,20])/2
  p[,4]<-(m2[,1]+m2[,17]+m2[,16])/3
  p[,5]<-(m2[,3]+m2[,9])/2
  p[,6]<-(m2[,6]+m2[,7]+m2[,4])/3
  p[,7]<-(m2[,2]+m2[,12])/2
  p[,8]<-m2[,5]
  p[,9]<-m2[,8]
  p[,10]<-m2[,15]
  x<-apply(p,2,mean)
  names(x)<-NULL
  D<-rep(0,10)
  for(j in 1:10){
    for(i in 1:L){
      D[j]<-D[j]+(p[i,j]-x[j])^2
    }
  }
  D<-(1/L)*D
  DD<-matrix(0,10,10)
  for(s in 1:10){
    for(t in 1:10){
      for(i in 1:(L-1)){
        DD[s,t]<-DD[s,t]+((p[i,s]-p[i+1,t])^2)/2
      }
    }
  }
  DD<-DD/(L-1)
  v<-c()
  for(i in 1:10){
    v<-c(v,DD[i,])
  }
  v<-c(v,D)
  v<-round(v,digits = 4)
  return(v)

}
