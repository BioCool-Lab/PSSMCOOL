#' single Average feature
#' @description This descriptor is a variant of the Average Block descriptor and is designed to group together rows
#'related to the same amino acid, thus considering domains of a sequence with similar conservation rates.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return feature vector of length 400
#' @references
#' L. Nanni, A. Lumini, and S. J. T. S. W. J. Brahnam, "An empirical study of different approaches for protein classification,"
#' vol. 2014, 2014.
#' @seealso \code{\link{Averag_Block}}
#' @export
#' @examples
#' w<-single_Average(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
single_Average<-function(pssm_name){
  x<-read.delim(pssm_name,skip = 2,sep = "",header = FALSE)
  x<-x[-1,-c(1,23:44)]
  d<-which(x=="Lambda")
  if(length(d)!=0){
    x<-x[-c(d:dim(x)[1]),]
  }
  colnames(x)<-NULL
  rownames(x)<-NULL
  P<-x[,1]
  P<-as.vector(P)
  E<-x[,-1]
  E<-as.matrix(E)
  mode(E)<-"integer"
  E<-1/(1+exp(-E))
  colnames(E)<-NULL
  v<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  SA<-rep(0,400)
  for(j in 1:20){
    az<-which(P==v[j])
    if (length(az)==0) {
      az<-1:(length(P))
    }
    for (i in 1:20) {
      SA[((j-1)*20+i)]<-mean(E[,i][az])
    }
  }
  SA<-round(SA,digits = 3)
  return(SA)
}
