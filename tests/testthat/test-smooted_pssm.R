context("checking that the length of smoothed_PSSM feature vector is equal to 4*220")

test_that("whether the smoothed_PSSM gives us the expected output",{
  ss<-smoothed_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),7,11,c(2,3,8,9))
  expect_equal(dim(ss),c(4,220))
})
