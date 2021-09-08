#' auto covariance transformation feature vector
#' @description The AC variable measures the correlation of the same property between two residues
#'separated by a distance of lg along the sequence
#' @param pssm_name name of the PSSM Matrix file
#' @param lg a parameter which indicates distance between two residues
#' @note in use of this function The lg parameter must be less than the length of the smallest sequence
#' in the database.
#' @import utils
#' @return feature vector which its length depends on parameter lg
#' @references
#' Dong, Q., Zhou, S. and Guan, J. (2009) A new taxonomy-based protein fold recognition approach
#' based on autocross-covariance transformation, Bioinformatics, 25, 2655-2662.
#' @export
#' @examples
#' q<-pssm_ac(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),17)
pssm_ac <- function(pssm_name,lg=18){ #lg smaler than shortest protein length in database
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
  AC<-matrix(0,nrow = lg,ncol = 20)
  g<-0
  for(t in 1:lg){
    for(j in 1:20){
      for (i in 1:(L-t)) {
        g<-g+(s[i,j]-sbar[j])*(s[i+t,j]-sbar[j])
      }
      AC[t,j]<-g/(L-t)
      g<-0
    }
  }
  vec<-c()
  for(i in 1:lg){
    vec<-c(vec,AC[i,])
  }
  vec<-round(vec,digits = 4)
  return(vec)
}
