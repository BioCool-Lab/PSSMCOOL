#' SOMA PSSM Feature
#' @description In this function each column can be viewed as a stochastic time series, and each PSSM contains 20
#'columns, in other words, each PSSM contains 20 stochastic time series and Second-order moving average (SOMA)
#'algorithm is applied to these columns to extract SOMA PSSM feature vector.
#' @param pssm_name name of PSSM file
#' @import utils
#' @return feature vector of length 160
#' @references
#' Y. Liang, S. J. J. o. M. G. Zhang, and Modelling, "Predict protein structural class by incorporating two different modes of evolutionary
#' information into Chou's general pseudo amino acid composition," vol. 78, pp. 110-117, 2017.
#' @export
#' @examples
#' w<-SOMA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
SOMA_PSSM<-function(pssm_name){
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
  somaPSSM<-rep(0,160)
  s<-0
  for(n in 2:9){
    for(j in 1:20){
      c2MA<-0
      Yhatn_i<-rep(0,L)
      for(i in n:L){
        Yhatn_i[i]<-sum(E[(i-(n-1)):i,j])
        c2MA<-c2MA+(E[i,j]-Yhatn_i[i])^2
      }
      c2MA<-1/(L-n)*c2MA
      s<-s+1
      somaPSSM[s]<-c2MA
    }
  }
  return(round(somaPSSM,digits = 4))
}
