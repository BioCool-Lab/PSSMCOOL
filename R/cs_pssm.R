#' CSP-SegPseP-SegACP feature vector
#' @description This feature vector is constructed  by fusing consensus sequence (CS), segmented PsePSSM, and
#'segmented auto-covariance transformation (ACT) based on PSSM. by consensus sequence a 40-dimensional feature
#'vector is obtained, in segmented PsePSSM group, by dividing PSSM Matrix to 2 and 3 segments a 380-dimensional
#'feature vector is obtained and in ACT group, similar to the previous group at first PSSM Matrix is divided
#'to 2 and 3 segments then a feature vector of length 280 is obtained.eventually by fusing these features a
#'700-dimensional feature vector is obtained.
#' @param pssm_name name of PSSM Matrix file
#' @param vec_name  a character that user imports to specify kind of feature vector which it can be varied between
#'four values
#' @import utils
#' @details If vec_name equals to "segmented_psepssm" then a feature vector of length 380 is obtained.
#'if vec_name equals to "segmented_acpssm" then a feature vector of length 280 is obtained, and if vec_name
#'equals to "cspssm" the obtained feature vector would be of length 40 eventually if vec_name equals to
#'"total" then feature vector would be of length 700.
#' @return feature vector that its length depends on the vec_name which user imports
#' @references
#' Y. Liang, S. Liu, S. J. C. Zhang, and m. m. i. medicine, "Prediction of protein structural
#' classes for low-similarity sequences based on consensus sequence and segmented PSSM," vol. 2015, 2015.
#' @export
#' @examples
#' A<-CS_PSe_PSSM(system.file("extdata", "C7GSI6.txt.pssm", package="PSSMCOOL"),"total")
CS_PSe_PSSM <- function(pssm_name,vec_name){
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
  names(L)<-NULL
  v<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  a<-rep(0,L)
  consensus<-c()
  for(i in 1:L){
    a[i]<-which.max(E[i,])
    consensus<-c(consensus,v[a[i]])
  }
  CSAAC<-c()
  CSCM<-c()
  i<-1
  while (i<=20) {
    t<-which(consensus==v[i])
    CSAAC[i]<-length(t)/L
    CSCM[i]<-sum(t)/(L*(L-1))
    i<-i+1
  }
  cspssm<-c(CSAAC,CSCM)
  n<-2
  L1<-round(L/n)
  L2<-L-L1
  a<-b<-matrix(0,5,20)
  o<-u<-v<-matrix(0,3,20)
  AC1<-AC2<-matrix(0,4,20)
  s<-ss<-0
  for(j in 1:20){
    a[1,j]<-sum(E[,j][1:L1])
    a[1,j]<-a[1,j]/L1
    b[1,j]<-sum(E[,j][(L1+1):L])
    b[1,j]<-b[1,j]/(L-L1)
    for(y in 1:4){
      for(i in 1:(L1-y)){
        s<-s+(E[i,j]-E[i+y,j])^2
      }
      a[y+1,j]<-s/(L1-y)
      s<-0
      for(i in (L1+1):(L-y)){
        s<-s+(E[i,j]-E[i+y,j])^2
      }
      b[y+1,j]<-s/(L-L1-y)
      s<-0
      for(i in 1:(L1-y)){
        ss<-ss+(E[i,j]-a[1,j])*(E[i+y,j]-a[1,j])
      }
      AC1[y,j]<-ss/(L1-y)
      ss<-0
      for(i in (L1+1):(L-y)){
        ss<-ss+(E[i,j]-b[1,j])*(E[i+y,j]-b[1,j])
      }
      AC2[y,j]<-ss/(L-L1-y)
    }
  }
  q1<-q2<-c()
  for(i in 1:4){
    q1<-c(q1,AC1[i,])
    q2<-c(q2,AC2[i,])
  }
  q3<-c(q1,q2)
  w1<-w2<-c()
  for(i in 1:5){
    w1<-c(w1,a[i,])
    w2<-c(w2,b[i,])
  }
  w3<-c(w1,w2)
  v1<-v2<-v3<-c()
  n<-3
  s<-ss<-0
  AC3<-AC4<-AC5<-matrix(0,2,20)
  L1<-round(L/3)
  L2<-2*L1
  L3<-(L-2*L1)
  for(j in 1:20){
    o[1,j]<-sum(E[,j][1:L1])
    o[1,j]<-o[1,j]/L1
    u[1,j]<-sum(E[,j][(L1+1):(2*L1)])
    u[1,j]<-u[1,j]/L1
    v[1,j]<-sum(E[,j][(2*L1 +1):L])
    v[1,j]<-v[1,j]/(L-2*L1)
    for(y in 1:2){
      for(i in 1:(L1-y)){
        s<-s+(E[i,j]-E[i+y,j])^2
      }
      o[y+1,j]<-s/(L1-y)
      s<-0
      for(i in (L1+1):(2*L1-y)){
        s<-s+(E[i,j]-E[i+y,j])^2
      }
      u[y+1,j]<-s/(L1-y)
      s<-0
      for(i in (2*L1 +1):(L-y)){
        s<-s+(E[i,j]-E[i+y,j])^2
      }
      v[y+1,j]<-s/(L-2*L1-y)
      for(i in 1:(L1-y)){
        ss<-ss+(E[i,j]-o[1,j])*(E[i+y,j]-o[1,j])
      }
      AC3[y,j]<-ss/(L1-y)
      ss<-0
      for(i in (L1+1):(2*L1-y)){
        ss<-ss+(E[i,j]-u[1,j])*(E[i+y,j]-u[1,j])
      }
      AC4[y,j]<-ss/(L1-y)
      ss<-0
      for(i in (2*L1+1):(L-y)){
        ss<-ss+(E[i,j]-v[1,j])*(E[i+y,j]-v[1,j])
      }
      AC5[y,j]<-ss/(L-2*L1-y)
      ss<-0
    }
  }
  q4<-q5<-q6<-c()
  for(i in 1:2){
    q4<-c(q4,AC3[i,])
    q5<-c(q5,AC4[i,])
    q6<-c(q6,AC5[i,])
  }
  q7<-c(q4,q5,q6)
  for(i in 1:3){
    v1<-c(v1,o[i,])
    v2<-c(v2,u[i,])
    v3<-c(v3,v[i,])
  }
  v4<-c(v1,v2,v3)
  segmented_psepssm<-c(w3,v4)
  segmented_acpssm<-c(q3,q7)
  total<-c(cspssm,segmented_psepssm,segmented_acpssm)
  vector_namse<-c("cspssm","segmented_psepssm","segmented_acpssm","total")
  switch(vec_name,vector_namse)
  dd<-eval(as.name(vec_name))
  dd<-round(dd,digits = 4)
  return(dd)
}
