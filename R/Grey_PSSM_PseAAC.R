#' grey pssm feature vector
#' @description This function produces a feature vector of length 100 which the first 20 components of this
#'vector is the normalized occurrence frequency of the native amino acids in the  protein. the next 20
#'components are mean of 20 PSSM columns and grey system model approach as elaborated in (Min et al. 2013)
#'is used to define the next 60 components.
#' @param pssm_name name of PSSM matrix file
#' @import utils
#' @return feature vector of length 100
#' @references
#' J.-L. Min, X. Xiao, and K.-C. J. B. r. i. Chou, "iEzy-Drug: A web server for identifying the interaction
#' between enzymes and drugs in cellular networking," vol. 2013, 2013.
#'
#' X. Xiao, M. Hui, and Z. J. T. J. o. m. b. Liu, "IAFP-Ense: an ensemble classifier for identifying
#' antifreeze protein by incorporating grey model and PSSM into PseAAC," vol. 249, no. 6, pp. 845-854, 2016.
#'
#' M. Kabir et al., "Improving prediction of extracellular matrix proteins using evolutionary information via
#' a grey system model and asymmetric under-sampling technique," vol. 174, pp. 22-32, 2018.
#' @export
#' @examples
#' as<-grey_pssm_pseAAC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
grey_pssm_pseAAC<-function(pssm_name){
  x<-read.delim(pssm_name,skip = 2,sep = "",header = FALSE)
  x<-x[-1,-c(1,23:44)]
  d<-which(x=="Lambda")
  if(length(d)!=0){
    x<-x[-c(d:dim(x)[1]),]
  }
  colnames(x)<-NULL
  rownames(x)<-NULL
  x<-as.matrix(x)
  m<-x
  protein_seq<-m[,1]
  protein_seq<-as.character(protein_seq)
  m2<-m[,-1]
  m2<-as.matrix(m2)
  mode(m2)<-"integer"
  m2<-1/(1+exp(-m2))
  L<-dim(m2)[1]
  feature_vect<-vector(mode = "numeric",length = 100)
  v<-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  for(i in 1:20){
    a<-which(protein_seq==v[i])
    feature_vect[i]<-length(a)/length(protein_seq)
    feature_vect[i]<-round(feature_vect[i],digits = 4)
  }
  ss<-apply(m2,2,mean)
  names(ss)<-NULL
  ss<-round(ss,digits = 4)
  feature_vect[21:40]<-ss
  B<-list()
  w<-list()
  s<-rep(0,60)
  for(j in 1:20){
    column3<-rep(1,L-1)
    column2<-rep(0,L-1)
    U<-rep(0,L-1)
    for(i in 1:(L-1)){
      column2[i]<- -(sum(m2[,j][1:i])+0.5*m2[i+1,j])
      U[i]<-m2[i+1,j]-m2[i,j]
    }
    U<-as.matrix(U)
    column1<- -m2[2:L,j]
    B[[j]]<-data.frame(column1,column2,column3)
    colnames(B[[j]])<-NULL
    rownames(B[[j]])<-NULL
    B[[j]]<-as.matrix(B[[j]])
    x<-t(B[[j]])%*%B[[j]]
    while(det(x)==0){
      diag(x)<-diag(x)+1
    }
    y<-t(B[[j]])%*%U
    x2<-solve(x)
    x2<-round(x2,digits = 4)
    w[[j]]<-x2%*%y
    s[3*j-2]<-feature_vect[j]*w[[j]][1,1]
    s[3*j-1]<-feature_vect[j]*w[[j]][2,1]
    s[3*j]<-feature_vect[j]*w[[j]][3,1]

  }
  feature_vect[41:100]<-s
  feature_vect<-round(feature_vect,digits = 4)
  return(feature_vect)
}
