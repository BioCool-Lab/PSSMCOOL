#' PSSM400 feature
#' @description This function firstly normalizes PSSM Matrix by formula:
#'\eqn{P-min(P)/max(P)-min(P)} then for any standard amino acid specifies its
#'position in protein sequence whereby a sub-matrix from PSSM corresponding
#'to these positions will be extracted, then for this sub-matrix computes
#'\code{\link[base]{colSums}}
#'of its columns to create a vector of length 20, eventually a feature vector of
#'length 400 will be obtained.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @note if a specific amino acid did not exist in protein then \code{\link[base]{colSums}} of whole
#'PSSM is computed.
#' @return feature vector of length 400
#' @export
#' @examples
#' q<-pssm400(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
pssm400<-function(pssm_name){
  x<-read.delim(pssm_name,skip = 2,sep = "",header = FALSE)
  x<-x[-1,-c(1,23:44)]
  d<-which(x=="Lambda")
  if(length(d)!=0){
    x<-x[-c(d:dim(x)[1]),]
  }
  colnames(x)<-NULL
  rownames(x)<-NULL
  p2<-x
  k2<-p2[,1]
  k2<-as.character(k2)
  p<-p2[,-1]
  p<-as.matrix(p)
  mode(p)<-"integer"
  p<-(p-min(p))/(max(p)-min(p))
  vec_final<-c()
  v<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  for(i in 1:20) {
    er<-which(k2==v[i])
    if (length(er)==0) {
      next
    } else if (length(er)==1) {
      v2<-p[k2==v[i],]
    }  else {
    b<-p[k2==v[i],]
    v2<-colSums(b)
    vec_final<-c(vec_final,v2)
    }
  }
  if (length(vec_final<400)) {
    vec_final<-c(vec_final,rep(0,(400-length(vec_final))))
  }
  vec_final<-round(vec_final,digits = 4)
  return(vec_final)
}
