context("checking that the length of SCSH feature vector is equal to 400 or 8000")

test_that("whether the SCSH gives us the expected output",{
  ax<-three_mer(3)
  ss<-scsh2(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),2)
  s1<-LETTERS[1:4]
  s2<-LETTERS[3:6]
  s<-data.frame(s1,s2)
  dc<-k_mers(s,3)
  expect_equal(length(ss),400)
  expect_equal(length(dc),16)
  expect_equal(dim(ax),c(2,8000))

})
