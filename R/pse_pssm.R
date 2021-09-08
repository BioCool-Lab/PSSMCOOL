#' pseudo position-specific scoring matrix feature
#' @description This feature vector is combination of \eqn{F_{PSSM}} feature vector and vector of
#'correlation factors correspond to 20 columns in PSSM Matrix. \eqn{F_{PSSM}} actually is mean of PSSM Matrix
#'columns of length 20.
#' @param pssm_name is the name of PSSM matrix file
#' @param g a parameter Which its size corresponds to the database used.
#' @import utils
#' @return feature vector of length 20+20*g
#' @references
#' D.-J. Yu et al., "Learning protein multi-view features in complex space," vol. 44, no. 5, pp. 1365-1379,
#' 2013.
#'
#' Chou, K.C. and Shen, H.B. (2007) MemType-2L: a web server for predicting membrane proteins and their
#' types by incorporating evolution information through Pse-PSSM, Biochemical and Biophysical Research
#' Communications, 360, 339-345.
#' @export
#' @examples
#' v<-pse_pssm(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
pse_pssm <- function(pssm_name,g=15){
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
  Ebar<-apply(E,2,mean)
  names(Ebar)<-NULL
  Ebar<-round(Ebar,digits = 4)
  G<-rep(0,20*g)
  for(lg in 1:g){
   for (j in 1:20) {
     k<-j+20*(lg-1)
    for(i in (1:(L-lg))){
      G[k]<-G[k]+(E[i,j]-E[i+lg,j])^2
    }
    G[k]<-(1/(L-lg))*G[k]
   }
  }

  v<-c(Ebar,G)
  v<-round(v,digits = 4)

  return(v)

}
