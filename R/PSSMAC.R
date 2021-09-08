#' PSSMAC feature
#' @description This feature, which stands for auto covariance transformation, for jth column calculates the average
#'of this column, and then subtracts the resulting number from the elements on the i and (i + g)th rows of this column,
#'and finally multiplies them. by changing the variable i from 1 to L-g, it calculates the sum of these, since the
#'variable j changes between 1 and 20, and the variable g between 1 and 10 eventually a feature vector of length 200
#'will be obtained.
#' @param pssm_name name of PSSM Matrix files
#' @import utils
#' @return feature vector of length 200
#' @references
#' L. Zou, C. Nan, and F. J. B. Hu, "Accurate prediction of bacterial type IV secreted effectors using amino acid
#' composition and PSSM profiles," vol. 29, no. 24, pp. 3135-3142, 2013.
#' @export
#' @examples
#' w<-PSSMAC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
PSSMAC <- function(pssm_name){
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
  s<-0
  x<-apply(p,2,mean)
  names(x)<-NULL
  x<-round(x,digits = 4)
  PSSM_AC<-matrix(0,20,10)
  for(g in 1:10){
    for(j in 1:20){
      for(i in 1:(L-g)){
        s<-s+(p[i,j]-x[j])*(p[i+g,j]-x[j])
      }
      PSSM_AC[j,g]<-s/(L-g)
      s<-0
    }
  }
  v<-c()
  for(j in 1:20){
    v<-c(v,PSSM_AC[j,])
  }
  v<-round(v,digits = 4)
  return(v)
}





