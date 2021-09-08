#' consunsus_sequence
#' @description This feature vector is constructed from PSSM Matrix as:
#'\eqn{\alpha(i)=argmax(P_{i,j})} where i varies between 1 and L and j between 1 and 20, L indicates protein length
#'and "arg" represents the argument of the maximum the ith base of the consensus sequence (CS) is then set to be the
#'\eqn{\alpha(i)}th amino acid in the amino acid alphabet and a consensus sequence is constructed.
#' @param pssm_name is the name of PSSM Matrix file
#' @import utils
#' @return consunsus sequence wich extracted from PSSM
#' @references
#' Y. Liang, S. Liu, S. J. C. Zhang, and m. m. i. medicine, "Prediction of protein structural
#' classes for low-similarity sequences based on consensus sequence and segmented PSSM," vol. 2015, 2015.
#' @export
#'
#' @examples
#' w<-consunsus_sequence(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
consunsus_sequence<-function(pssm_name){
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
  L<-dim(m)[1]
  aaVect<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  d<-c()
  for(i in 1:L){
    d[i]<-which.max(m[i,])
  }
  consunsus<-aaVect[d]
  consunsus<-paste(consunsus,collapse = "")
  return(consunsus)
}

