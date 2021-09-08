#' Disulfide connectivity feature
#' @description This feature is used to predict the disulfide bond within a protein.
#' @param pssm_name name of PSSM Matrix file
#' @import utils
#' @details  For the purpose of predicting disulfide bond in protein at first, the total
#'number of cysteine amino acids in the protein sequence is counted and their position in the protein sequence is
#'identified. Then, using a sliding window with length of 13, moved on the PSSM matrix from top to bottom so that
#'the middle of the window is on the amino acid cysteine, then the rows below the matrix obtained from the PSSM
#'matrix with dimension of 13 x 20 are placed next to each other to get a feature vector with a length of
#'260 = 20 * 13 per cysteine, and if the position of the first and last cysteine in the protein sequence is such
#'that the middle of sliding window is not on cysteine residue when moving on PSSM Matrix, then the required
#'number of zero rows from top and bottom is added to the PSSM matrix to achieve this goal.Thus, for every
#'cystine amino acid presented in protein sequence, a feature vector with a length of 260 is formed.Then all
#'the pairwise combinations of these cysteines is wrote in the first column of a table, and in front of each of these
#'pairwise combinations, the corresponding feature vectors are sticked together to get a feature vector of length 520
#'for each of these compounds.Finally, the table obtained in this way will have the number of rows equal to the
#'number of all pairwise combinations of these cysteines and the number of columns will be equal to 521
#'(the first column includes the name of these pair combinations). And it is easy to divide this table into
#'training and testing data and predict the desired disulfide bonds between cysteines.
#' @return a table with number of all cyctein pairs in rows and 521 columns correspond to feature vector length.
#' @references
#' D.-J. Yu et al., "Disulfide connectivity prediction based on modelled protein 3D structural information and
#' random forest regression," vol. 12, no. 3, pp. 611-621, 2014.
#'
#' N. J. Mapes Jr, C. Rodriguez, P. Chowriappa, S. J. C. Dua, and s. b. journal, "Residue adjacency matrix based
#' feature engineering for predicting cysteine reactivity in proteins," vol. 17, pp. 90-100, 2019.
#' @export
#' @examples
#' aq<-disulfid(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
disulfid<-function(pssm_name){
  x<-read.delim(pssm_name,skip = 2,sep = "",header = FALSE)
  x<-x[-1,-c(1,23:44)]
  d<-which(x=="Lambda")
  if(length(d)!=0){
    x<-x[-c(d:dim(x)[1]),]
  }
  colnames(x)<-NULL
  rownames(x)<-NULL
  x<-as.matrix(x)
  k<-x[,1]
  k<-as.character(k)
  p<-x[,-1]
  mode(p)<-"integer"
  p<-1/(1+exp(-p))
  L<-dim(p)[1]
  v<-which(k=="C")
  if(length(v)==0|length(v)==1){
    stop("there is no disulfid bond")
  } else {
    if(v[1]<7){
      h<-matrix(0,(7-v[1]),20)
      p<-rbind(h,p)
      v<-v+(7-v[1])
    }
    if(v[length(v)]>(L-6)){
      tt<-matrix(0,(v[length(v)]-(L-6)),20)
      p<-rbind(p,tt)
    }
    w<-matrix(0,length(v),260)
    for(i in v){
      x<-p[(i-6):(i+6),]
      a<-which(v==i)
      ss<-c()
      for(j in 1:dim(x)[1]){
        ss<-c(ss,x[j,])
      }
      w[a,]<-ss
    }
    colnames(w)<-NULL
    az<-c()
    aw<-c()
    for(i in 1:(dim(w)[1]-1)){
      for(j in (i+1):dim(w)[1]){
        az<-c(az,paste0("c",i,"c",j))
        aw<-rbind(aw,c(w[i,],w[j,]))
      }
    }
  }
  aw<-round(aw,digits = 4)
  aw<-as.data.frame(aw)
  aw<-cbind(az,aw)
  colnames(aw)<-1:521
  return(aw)
}
