#' D-FPSSM and SF-PSSM feature vectors
#' @description This function produces list of two feature vectors named D-FPSSM and S-FPSSM
#'which then used by FPSSM2 function to construct feature vector of length 100 for each pair of
#'proteins which then used for protein-protein interaction prediction in each dataset.
#' @param pssm_name name of PSSM Matrix file
#' @param hk a parameter that indicates which amino acid alphabet must be used
#' @import utils
#' @return two feature vectors of different length which is used in later steps.
#' @references
#' Zahiri, J., et al. (2013) PPIevo: protein-protein interaction prediction from PSSM based
#' evolutionary information, Genomics, 102, 237-242.
#' @export
#' @examples
#' q<-FPSSM(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),8)
#'
FPSSM <- function(pssm_name,hk){                              #### hk must be 8 or 20
  x<-read.delim(pssm_name,skip = 2,sep = "",header = FALSE)
  x<-x[-1,-c(1,23:44)]
  d<-which(x=="Lambda")
  if(length(d)!=0){
    x<-x[-c(d:dim(x)[1]),]
  }
  colnames(x)<-NULL
  rownames(x)<-NULL
  x<-as.matrix(x)
  k2<-x[,1]
  k2<-as.character(k2)
  p<-x[,-1]
  mode(p)<-"integer"
  p[is.na(p)]<-0
  p[p<0]<-0
  if(hk==8){
    Group<-list(Grp1<-c("A","E"),Grp2<-c("R","Q","K","H"),Grp3<-c("N","D","S","T")
                ,Grp4<-c("G"),Grp5<-c("P"),Grp6<-c("I","L","M","F","V"),Grp7<-c("W","Y"),Grp8<-"C")
    p<-p[,c(1,7,2,6,12,9,3,4,16,17,8,15,10,11,13,14,20,18,19,5)]
    s<-matrix(0,8,20)
    D<-apply(p,2,sum)
    names(D)<-NULL
    h1<-c(sum(D[1:2]),sum(D[3:6]),sum(D[7:10]),D[11],D[12],sum(D[13:17]),sum(D[18:19]),D[20])
    D<-h1
    q<-c()
    for(i in 1:8){
      for(j in 1:20){
        s[i,j]<-sum(p[,j][k2 %in% Group[[i]]])
      }
      h<-c(sum(s[i,][1:2]),sum(s[i,][3:6]),sum(s[i,][7:10]),s[i,11],s[i,12],sum(s[i,][13:17]),sum(s[i,][18:19]),s[i,20])
      q<-c(q,h)
    }
    S<-q
    }
  else if(hk==20){
    S<-c()
    s<-matrix(0,20,20)
    D<-apply(p,2,sum)
    names(D)<-NULL
    D<-round(D,digits =4)
    v<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
    for(i in 1:20){
      for(j in 1:20){
        s[i,j]<-sum(p[,j][k2==v[i]])
      }
      S<-c(S,s[i,])
    }
  }
  else stop("ERROR: entered Number must be 8 or 20")

  S<-round(S,digits =3)
  return(list(D,S))
}
#######################################################################################################
#' Mixture of Two FPSSM Features
#' @description This function takes two PSSM files as argument and uses FPSSM function for making feature
#'vector of length 100 correspond to this pair of proteins.
#' @param pssm_name1 The name of first PSSM Matrix file
#' @param pssm_name2 The name of second PSSM Matrix file
#' @param hk a parameter that indicates which amino acid alphabet must be used
#' @importFrom stats var
#' @importFrom infotheo mutinformation
#' @importFrom utils read.delim
#' @return Feature vector of length 100
#' @references
#' Zahiri, J., et al. (2013) PPIevo: protein-protein interaction prediction from PSSM based
#' evolutionary information, Genomics, 102, 237-242.
#' @seealso \code{\link[entropy]{entropy}}
#'
#' \code{\link[infotheo]{mutinformation}}
#' @export
#' @examples
#'  s1<-system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL")
#'  s2<-system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL")
#'  s<-FPSSM2(s1,s2,8)
FPSSM2<-function(pssm_name1,pssm_name2,hk){
  requireNamespace("infotheo",quietly = TRUE)
  if(hk==20){
  f100<-c()
  L1<-FPSSM(pssm_name1,hk)
  S1<-L1[[2]]
  D1<-L1[[1]]
  L2<-FPSSM(pssm_name2,hk)
  S2<-L2[[2]]
  D2<-L2[[1]]
  S3<-abs(S1-S2)
  v1=v2=v3<-rep(0,20)
  for(i in 1:20){
    v1[i]<-min(S3[((i-1)*20+1):(i*20)])
    v2[i]<-var(S3[((i-1)*20+1):(i*20)])
    v3[i]<-infotheo::mutinformation(S1[((i-1)*20+1):(i*20)],S2[((i-1)*20+1):(i*20)])
    f100<-c(f100,v1[i],v2[i],v3[i])
  }
  f100<-c(f100,D1,D2)
  f100<-round(f100,digits =2)
  return(f100)
  }
  else if(hk==8){
    f40<-c()
    L1<-FPSSM(pssm_name1,hk)
    S1<-L1[[2]]
    D1<-L1[[1]]
    L2<-FPSSM(pssm_name2,hk)
    S2<-L2[[2]]
    D2<-L2[[1]]
    S3<-abs(S1-S2)
    v1=v2=v3<-rep(0,8)
    for(i in 1:8){
      v1[i]<-min(S3[((i-1)*8+1):(i*8)])
      v2[i]<-var(S3[((i-1)*8+1):(i*8)])
      v3[i]<-infotheo::mutinformation(S1[((i-1)*8+1):(i*8)],S2[((i-1)*8+1):(i*8)])
      f40<-c(f40,v1[i],v2[i],v3[i])
    }
    f40<-c(f40,D1,D2)
    f40<-round(f40,digits =2)
    return(f40)
  }
  else{
    stop("ERROR: parameter value must be 8 or 20")
  }
}
