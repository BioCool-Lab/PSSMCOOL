#' MBMGACPSSM feature
#' @description In this function three different autocorrelation descriptors based on PSSM are adopted,
#'which include: normalized Moreau-Broto autocorrelation, Moran autocorrelation and Geary autocorrelation
#'descriptors.Autocorrelation descriptor is a powerful statistical tool and defined based on the distribution
#'of amino acid properties along the sequence, which measures the correlation between two residues separated
#'by a distance of d in terms of their evolution scores.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @return feature vector of length 560
#' @references
#' Y. Liang, S. Liu, S. J. M. C. i. M. Zhang, and i. C. Chemistry, "Prediction of protein structural class
#' based on different autocorrelation descriptors of position-specific scoring matrix,"
#' vol. 73, no. 3, pp. 765-784, 2015.
#' @export
#' @examples
#' w<-MBMGACPSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
MBMGACPSSM<-function(pssm_name){
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
  P_Bar<-apply(p,2,mean)
  names(P_Bar)<-NULL
  L<-dim(p)[1]
  s<-0
  P_MBACP<-matrix(0,9,20)
  for(d in 1:9){
    for(j in 1:20){
      for(i in 1:(L-d)){
        s<-s+p[i,j]*p[i+d,j]
      }
      a<-s/(L-d)
      a<-round(a,digits = 3)
      P_MBACP[d,j]<-a
      s<-0
    }
  }
  v1<-c()
  for(i in 1:9){
    v1<-c(v1,P_MBACP[i,])
  }
  s<-0
  s2<-0
  P_MACP<-matrix(0,9,20)
  for(d in 1:9){
    for(j in 1:20){
      for(i in 1:(L-d)){
        s<-s+(p[i,j]-P_Bar[j])*(p[i+d,j]-P_Bar[j])
      }
      s<-1/(L-d)*s
      for(t in 1:L){
        s2<-s2+(p[t,j]-P_Bar[j])^2
      }
      s2<-1/L*s2
      ee<-s/s2
      P_MACP[d,j]<-round(ee,digits = 3)
      s=s2<-0
    }
  }
  v2<-c()
  for(i in 1:9){
    v2<-c(v2,P_MACP[i,])
  }
  s=s2=0
  P_GACP<-matrix(0,9,20)
    for(d in 1:9){
      for(j in 1:20){
        for(i in 1:(L-d)){
          s<-s+(p[i,j]-p[i+d,j])^2
        }
        s<-1/2*(L-d)*s
        for(t in 1:L){
          s2<-s2+(p[t,j]-P_Bar[j])^2
        }
        s2<-1/(L-1)*s2
        aa<-s/s2
        P_GACP[d,j]<-round(aa,digits = 3)
        s=s2<-0
      }
    }
  v3<-c()
  v3<-c(v3,P_GACP)
  MBMGAC_PSSM<-c(P_Bar,v1,v2,v3)
  MBMGAC_PSSM<-round(MBMGAC_PSSM,digits = 3)
  return(MBMGAC_PSSM)

}
