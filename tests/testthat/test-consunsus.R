context("checking that the length of consunsus_sequence is equal to number of PSSM Matrix rows")

test_that("whether the consunsus_sequence function gives us the expected output",{
  ss<-consunsus_sequence(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  x<-read.delim(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),skip = 2,sep = "",header = FALSE)
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
  expect_equal(nchar(ss),dim(x)[1])
})
