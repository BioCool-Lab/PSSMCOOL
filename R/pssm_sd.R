#' PSSM-SD feature
#' @description In this feature, by considering a specific column, at first sum of all components in this column
#'is denoted by "L", then starting from the first row in this column, the components are added together to
#'reaching a value less than or equal to 25% of L, and the number of components to get such a sum in this column
#'is calculated and stored . In the next step, the same work is done starting from the first row to reaching a
#'value less than or equal to 50% of L, and the number of components is saved. In the next step, the same work
#'is done started from the last row, To reaching 25% and 50% of L, and number of components is saved as before.
#'By appending these saved numbers together for each column, a vector of length 4 is obtained. If this is done
#'for all the columns and the obtained vectors are connected to each other, for each protein, a feature vector
#'of length 80 is obtained which its name is PSSM-SD.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return feature vector of length 80
#' @references
#' A. Dehzangi, K. Paliwal, J. Lyons, A. Sharma, A. J. I. A. T. o. C. B. Sattar, and Bioinformatics,
#' "A segmentation-based method to extract structural and evolutionary features for protein fold
#' recognition," vol. 11, no. 3, pp. 510-519, 2014.
#' @export
#' @examples
#' ww<-PSSM_SD(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
PSSM_SD <- function(pssm_name){
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
  #p<-1/(1+exp(-p))
  L<-dim(p)[1]
  t<-apply(p,2,sum)
  names(t)<-NULL
  t<-round(t,digits = 4)
  I<-matrix(0,4,20)
  for(j in 1:20){
    i<-1;r<-L
    s<-p[i,j]
    sr<-p[r,j]
    while(s>(1/4)*t[j]) {
      i<-i+1
      s<-s+p[i,j]
    }
    I[1,j]<- i
    i<-1
    s<-p[i,j]
    while (s>(1/2)*t[j]) {
      i<-i+1
      s<-s+p[i,j]
    }
    I[2,j]<- i
    while (sr>(1/4)*t[j]) {
      r<-r-1
      sr<-sr+p[r,j]
    }
    I[3,j]<- r
    r<-L
    sr<-p[r,j]
    while (sr>(1/2)*t[j]) {
      r<-r-1
      sr<-sr+p[r,j]
    }
    I[4,j]<- r
  }
  vv<-c()
  I<-round(I,digits = 4)
  for(h in 1:4){
    vv<-c(vv,I[h,])
  }
  return(list(I,vv))
}


#############################################################################################################
#' PSSM-Seg feature vector
#' @description This feature vector uses PSSM-SD to produce Segmented Auto Covariance Features.
#' @param pssm_name name of PSSM Matrix file
#' @param m a parameter between 1 and 11
#' @import utils
#' @return feature vector of length 100
#' @references
#' A. Dehzangi, K. Paliwal, J. Lyons, A. Sharma, A. J. I. A. T. o. C. B. Sattar, and Bioinformatics,
#' "A segmentation-based method to extract structural and evolutionary features for protein fold
#' recognition," vol. 11, no. 3, pp. 510-519, 2014.
#' @seealso \code{\link{PSSM_SD}}
#' @export
#' @examples
#' q<-pssm_seg(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),3)
pssm_seg <- function(pssm_name,m=4){
  I<-PSSM_SD(pssm_name)[[1]]
  if(min(I)< m){
    m<-(min(I)-1)
  }
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
  mt<-apply(p,2,mean)
  names(mt)<-NULL
  mt<-round(mt,digits = 4)
  seg<-matrix(0,4,20)
  s<-0
  for(t in 1:4){
    for(j in 1:20){
      for(i in 1:(I[t,j]-m)){
        s<-s+(p[i,j]-mt[j])*(p[i+m,j]-mt[j])
      }
      seg[t,j]<-(1/(I[t,j]-m))*s
      s<-0
    }
  }
  s<-0
  pssm_ac<-rep(0,20)
  for(j in 1:20){
    for(i in 1:(L-m)){
      s<-s+(p[i,j]-mt[j])*(p[i+m,j]-mt[j])
    }
    pssm_ac[j]<- (1/(L-m))*s
    s<-0
  }
  pseg<-c()
  for(j in 1:20){
    pseg<-c(pseg,seg[,j])
  }
  PSSM_SAC<-c(pseg,pssm_ac)
  PSSM_SAC<-round(PSSM_SAC,digits = 4)
  return(PSSM_SAC)
}

