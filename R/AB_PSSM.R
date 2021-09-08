#' AB-PSSM and RPM-PSSM feature vector
#' @description This feature consists of two types of feature vectors. at first, each protein sequence is divided into
#'20 equal parts, each of which is called a block, and in each block the row vectors of the PSSM matrix related
#'to that block are added together and The resulting final vector is divided by the length of that block, which is
#'5% of the total length of the protein.
#'Finally, by placing these 20 vectors side by side, the first feature vector of length 400 is obtained. The second feature
#'for each amino acid in each column is the average of the positive numbers in that column and for each block, and these 20
#'numbers, corresponding to 20 blocks, are placed next to each other, and therefore For each of the 20 types of amino acids,
#'a vector of length 20 is obtained, and by placing these together, the vector of the second feature, length 400, is obtained.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return two feature vectors AB-PSSM and RPM-PSSM each of length 400
#' @references
#' Jeong, J.C., Lin, X. and Chen, X.W. (2011) On position-specific scoring matrix for protein function prediction
#' , IEEE/ACM transactions on computational biology and bioinformatics / IEEE, ACM, 8, 308-315.
#' @export
#' @examples
#' zz<- AB_PSSM(system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL"))
AB_PSSM <- function(pssm_name){
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
  s<-floor(L/20)
  sm<-rep(0,20)
  p<-0;n<-0
  dd<-c()
  f<-g<-matrix(0,20,20) # g maked for feature set 2
  for(j in 1:20){
    for(t in 1:19){
      r1<-(1+(t-1)*s)
      r2<-(s+(t-1)*s)
      for(i in r1:r2){
        sm<-sm+m2[i,]
        if(m2[i,j]>0){
          p<-p+m2[i,j]
          n<-n+1
        }
      }
      if(p==0){
        dd[t]<-0
      }
      else{
        dd[t]<-p/n
      }
      p<-0;n<-0
      f[,t]<-sm/s
      sm<-rep(0,20)

    }
    dd[20]<-0
    g[,j]<-dd
  }
  e<-(19*s+1)
  for(i in e:L){
    sm<-sm+m2[i,]
  }
  f[,20]<-sm/(L-19*s)
  for(j in 1:20){
    for(i in e:L){
      if(m2[i,j]>0){
        p<-p+m2[i,j]
        n<-n+1

      }
    }
    if(p==0){
      g[20,j]<-0
    }
    else{
      g[20,j]<-p/n
    }
    p<-n<-0
  }
  v<-w<-c()
  for(i in 1:20){
    v<-c(v,f[,i])
    w<-c(w,g[,i])
  }
  v<-round(v,digits = 4)
  w<-round(w,digits = 4)
 return(list(v,w))
}




