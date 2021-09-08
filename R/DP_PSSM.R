#' DP_PSSM feature vector
#' @description This feature results from the connection of two vectors. The vector is the first feature of a
#'vector with a length of 40, which calculates the average of positive and negative values for each column
#'separately and puts them together. in the second feature vector, correspond to each column the difference
#'between the numbers in the rows that have distance of k is calculated, and then the square average for the
#'differences that are positive is calculated, and the same action  for the differences that are negative is
#'performed. since k varies between 1 and \eqn{\alpha}, and because the value of \eqn{\alpha} in this function is equal to 2,
#'the length of the second feature vector will be 80, which by merging with the first feature vector, the total
#'feature vector of length 120 will be obtained.
#' @param pssm_name name of PSSM matrix file
#' @param a fixed parameter that user chooses which usually equals to 2
#' @importFrom stats sd
#' @import utils
#' @return feature vector of length 120
#' @references
#' Juan, E.Y., et al. (2009) Predicting Protein Subcellular Localizations for Gram-Negative Bacteria using DP-PSSM
#' and Support Vector Machines. Complex, Intelligent and Software Intensive Systems, 2009. CISIS'09. International
#' Conference on. IEEE, pp. 836-841.
#' @export
#' @examples
#' ss<-DP_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
DP_PSSM <- function(pssm_name,a=2){
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
  m2<-x
  m2<-1/(1+exp(-m2))
  L<-dim(m2)[1]
  mt<-matrix(0,L,20)
  for(i in 1:L){
    if(sd(m2[i,])==0){
      mt[i,]<-m2[i,]
    }else{
      mt[i,]<-(m2[i,]-mean(m2[i,]))/sd(m2[i,])
    }
  }
  t1<-1:40
  for(j in 1:20){
    p<-mt[,j][mt[,j]>0]
    n<-mt[,j][mt[,j]<0]
    t1[2*j-1]<-mean(p)
    t1[2*j]<-mean(n)

  }
  pp<-matrix(0,L-a,20)
  np<-c() #number of positive differences in column j
  nn<-c() #number of negative differences in column j

  for(j in 1:20){
    for(i in 1:(L-a)){
      pp[i,j]<-mt[i,j]-mt[i+a,j]
    }
    ss<-pp[,j][pp[,j]>0]
    ff<-pp[,j][pp[,j]<0]
    np[j]<-length(ss)
    nn[j]<-length(ff)
  }
  x<-c()
  smp<-0
  smn<-0
  deljp<-rep(0,a)
  deljn<-rep(0,a)
  G<-c()
  for(j in 1:20){
    for(k in 1:a){
      for(i in 1:(L-k)){
        if((mt[i,j]-mt[i+k,j])>0){
          smp<-smp+(mt[i,j]-mt[i+k,j])^2
        }
        if((mt[i,j]-mt[i+k,j])<0){
          smn<-smn+(mt[i,j]-mt[i+k,j])^2
        }
      }
      deljp[k]<-smp/np[j]
      x<-c(x,deljp[k])
      deljn[k]<-(-smn/nn[j])
      x<-c(x,deljn[k])
      smp<-0
      smn<-0
    }
    G<-c(x,G)
    x<-c()
  }
  H<-c(t1,G)
  H<-round(H,digits = 4)
  return(H)
}


