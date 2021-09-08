#' trigrame feature vector
#' @description This feature vector is 8000-dimentional feature vector wich is computed from  tri-gram probability matrix
#'T obtained from PSSM Matrix.to achieve this purpose elements in three successive rows and arbitrary columns  are multiplied
#'together then these results are added together by changing variable i from 1 to L-1, which i is counter of row and L
#'indicates protein length. since there are 20 columns thus final feature vector would be of length 8000.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return feature vector of lenght 8000
#' @references
#' Paliwal, K.K., et al. (2014) A tri-gram based feature extraction technique using linear probabilities of position
#' specific scoring matrix for protein fold recognition, IEEE transactions on nanobioscience, 13, 44-50
#' @export
#' @examples
#' as<-trigrame_pssm(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GSI6.txt.pssm"))
trigrame_pssm<-function(pssm_name){
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
  t<-array(0,dim = c(20,20,20))
  k<-1
  vec<-vector(mode = "numeric",length = 8000)
  for(m in 1:20){
    for(n in 1:20){
      for(r in 1:20){
        for(i in 1:(L-2)){
          t[m,n,r]<-t[m,n,r]+p[i,m]*p[i+1,n]*p[i+2,r]
        }
        vec[k]<-t[m,n,r]
        k<-k+1
      }
    }
  }
  return(round(vec,digits = 4))
}

#system.time({trigrame_pssm(paste0(system.file("extdata",package="PSSMFeatures"),"/C7GSI6.txt.pssm"))})
