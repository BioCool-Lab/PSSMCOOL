#' kiderafactor feature
#' @description For product of this feature vector similar to smoothed_PSSM feature, firstly PSSM Matrix is
#'smoothed by appending zero vectors to its head and tail and utilizing sliding window of size odd, then this
#'smoothed PSSM Matrix is condensed by the Kidera factors to produce feature vector for each residue.
#' @param pssm_name name of PSSM Matrix file
#' @param v vector of amino acids positions which we want to produce feature vector for them.
#' @import utils
#' @return matrix of feature vectors
#' @references
#' C. Fang, T. Noguchi, H. J. I. j. o. d. m. Yamana, and bioinformatics, "Condensing position-specific scoring
#' matrixs by the Kidera factors for ligand-binding site prediction," vol. 12, no. 1, pp. 70-84, 2015.
#' @seealso \code{\link{smoothed_PSSM}}
#' @export
#' @examples
#' w<-kiderafactor(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),c(2,3,8,9))
kiderafactor<-function(pssm_name,v=NULL){
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
  m<-1/(1+exp(-m))
  L<-dim(m)[1]
  smoothed_PSSM<-matrix(0,L,20)
  h<-matrix(0,3,20)
  k<-matrix(0,3,20)
  m<-rbind(h,m,k)
  for(i in 1:L){
    E<-data.frame(m[i,],m[i+1,],m[i+2,],m[i+3,],m[i+4,],m[i+5,],m[i+6,])
    smoothed_PSSM[i,]<-rowSums(E)/7
  }
  v1<-c(-1.56,0.22,1.14,0.58,0.12,-0.47,-1.45,1.46,-0.41,-0.73,-1.04,-0.34,-1.40,-0.21,2.06,0.81,0.26,0.30,1.38,-0.74)
  v2<-c(-1.67,1.27,-0.07,-0.22,-0.89,0.24,0.19,-1.96,0.52,-0.16,0.00,0.82,0.18,0.98,-0.33,-1.08,-0.70,2.10,1.48,-0.71)
  v3<-c(-0.97,1.37,-0.12,-1.58,0.45,0.07,-1.61,-0.23,-0.28,1.79,-0.24,-0.23,-0.42,-0.36,-1.15,0.16,1.21,-0.72,0.80,2.04)
  v4<-c(-0.27,1.87,0.81,0.81,-1.05,1.10,1.17,-0.16,0.28,-0.77,-1.10,1.70,-0.73,-1.43,-0.75,0.42,0.63,-1.57,-0.56,-0.40)
  v5<-c(-0.93,-1.70,0.18,-0.92,-0.71,1.10,-1.31,0.10,1.61,-0.54,-0.55,1.54,2.00,0.22,0.88,-0.21,-0.10,-1.16,-0.00,0.50)
  v6<-c(-0.78,0.46,0.37,0.15,2.41,0.59,0.40,-0.11,1.01,0.03,-2.05,-1.62,1.52,-0.81,-0.45,-0.43,0.21,0.57,-0.68,-0.81)
  v7<-c(-0.20,0.92,-0.09,-1.52,1.52,0.84,0.04,1.32,-1.85,-0.83,0.96,1.15,0.26,0.67,0.30,-1.89,0.24,-0.48,-0.31,-1.07)
  v8<-c(-0.08,-0.39,1.23,0.47,-0.69,-0.71,0.38,2.36,0.47,0.51,-0.76,-0.08,0.11,1.10,-2.30,-1.15,-1.15,-0.40,1.03,0.06)
  v9<-c(0.21,0.23,1.10,0.76,1.13,-0.03,-0.35,-1.66,1.13,0.66,0.45,-0.48,-1.27,1.71,0.74,-0.97,-0.56,-2.30,-0.05,-0.46)
  v10<-c(-0.48,0.93,-1.73,0.70,1.10,-2.33,-0.12,0.46,1.63,-1.78,0.93,0.60,0.27,-0.44,-0.28,-0.23,0.19,-0.60,0.53,0.65)
  kidera_table<-data.frame(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10)
  kidera_table<-as.matrix(kidera_table)
  colnames(kidera_table)<-NULL
  mh<-matrix(0,8,20)
  mk<-matrix(0,8,20)
  condenced_pssm<-rbind(mh,smoothed_PSSM,mk)
  d<-9
  a<-1
  s<-0
  x<-matrix(0,L,170)
  w1=w2<-c()
  for(i in (d-8):(d+8)){
    for(p in 1:10){
      for(j in 1:20){
        s<-s+condenced_pssm[i,j]*kidera_table[j,p]
      }
      w2[p]<-s
      s<-0
    }
    w1<-c(w1,w2)
    w2<-c()
  }
  x[a,]<-w1
  d<-i+1
  a<-a+1
  while(d<=(L+16)){
    i==d
    s<-0
    for(p in 1:10){
      for(j in 1:20){
        s<-s+condenced_pssm[i,j]*kidera_table[j,p]
      }
      w2[p]<-s
      s<-0
    }
    w1<-c(w1[11:170],w2)
    x[a,]<-w1
    a<-a+1
    w2<-c()
    d<-d+1
  }
  if(length(v)!=0){
    y<-x[v,]
    y<-as.data.frame(y)
    rownames(y)<-v
    colnames(y)<-1:170
  }
  else{
    y<-x
    y<-as.data.frame(y)
    rownames(y)<-1:L
    colnames(y)<-1:170
  }
  return(round(y,digits = 3))
}






