#' DMACA-PSSM feature
#' @description In this feature each column of PSSM Matrix, can be regarded as a time series. Each PSSM
#'contains 20 columns Hence, each PSSM can be considered as 20 time series.The detrended moving-average
#'cross-correlation analysis (DMCA) is developed to measure the level of cross-correlation between two
#'non-stationary time series by fusing the detrended cross-correlation analysis (DCCA) and the detrended
#'moving average(DMA).this function utilises this algorithm for each column and each pair of columns to produce
#'a feature vector of length 290.
#' @param pssm_name name of PSSM Matrix file
#' @param n A parameter called the window size that must be smaller than the length of the sequence
#' @note parameter n must be equal or greater than 3 and equal or less then L which L is length of protein
#' @import utils
#'
#' @return feature vector of length 210
#' @references
#' Y. Liang, S. Zhang, S. J. S. Ding, and Q. i. E. Research, "Accurate prediction of Gram-negative bacterial secreted protein
#' types by fusing multiple statistical features from PSI-BLAST profile," vol. 29, no. 6, pp. 469-481, 2018.
#'
#' Y. Liang and S. J. A. b. Zhang, "Prediction of apoptosis protein’s subcellular localization by fusing two different
#' descriptors based on evolutionary information," vol. 66, no. 1, pp. 61-78, 2018.
#' @export
#' @examples
#' as<-DFMCA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7)
DFMCA_PSSM<-function(pssm_name,n){
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
  DMAA_n<-rep(0,20)
  DMCA_vec<-rep(0,190)
  s<-0
  for(i in 1:20){
    X<-rep(0,L)
    for(k in 1:L){
      X[k]<-sum(E[1:k,i])
    }
    Xhatk_n<-rep(0,(L-(n-1)))
    for(t in 1:(L-(n-1))){
      Xhatk_n[t]<-(1/n)*sum(X[t:(t+(n-1))])
    }
    for(j in i:20){
      Y<-rep(0,L)
      for(k in 1:L){
        Y[k]<-sum(E[1:k,j])
      }
      DMCA_n<-0
      Yhatk_n<-rep(0,(L-(n-1)))
      for(t in 1:(L-(n-1))){
        Yhatk_n[t]<-(1/n)*sum(Y[t:(t+(n-1))])
        DMCA_n<-DMCA_n+(X[t]-Xhatk_n[t])*(Y[t]-Yhatk_n[t])
      }
      DMCA_n<-1/(L-(n-1))*DMCA_n
      if(j==i){
        DMAA_n[i]<-DMCA_n
      }
      else{
        s<-s+1
        DMCA_vec[s]<-DMCA_n
      }
    }
  }

  DMAA_n<-round(DMAA_n,digits = 4)
  DMCA_vec<-round(DMCA_vec,digits = 4)
  return(c(DMAA_n,DMCA_vec))
}

