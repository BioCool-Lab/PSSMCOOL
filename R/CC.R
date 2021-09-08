#' Cross covarianse feature vector
#' @description The CC variable measures the correlation of two different properties between two residues
#'separated by a distance of t along the sequence.
#' @param pssm_name name of PSSM Matrix file
#' @param t shortest protein length in dataset minus one
#' @import utils
#' @return feature vector of length 380
#' @references
#' Dong, Q., Zhou, S. and Guan, J. (2009) A new taxonomy-based protein fold recognition approach
#' based on autocross-covariance transformation, Bioinformatics, 25, 2655-2662.
#' @export
#' @examples
#' aa<-pssm_cc(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),18)
pssm_cc <- function(pssm_name,t){
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
  s<-x
  s<-1/(1+exp(-s))
  L<-dim(s)[1]
  sbar<-apply(s,2,mean)
  names(sbar)<-NULL
  sbar<-round(sbar,digits = 4)
  h<-1:20
  g<-0
  n<-1
  v<-vector(mode = "numeric",length = 380)
    for(j1 in h){
      for(j2 in h[-j1]){
        for(i in 1:(L-t)){
          g<-g+(s[i,j1]-sbar[j1])*(s[i+t,j2]-sbar[j2])
        }
        v[n]<-g/(L-t)
        g<-0
        n<-n+1
      }

    }
  v<-round(v,digits = 4)
  return(v)
}
