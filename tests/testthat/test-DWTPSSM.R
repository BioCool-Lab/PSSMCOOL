context("checking that the length of dwt_PSSM feature vector is equal to 80")

test_that("whether the dwt_PSSM gives us the expected output",{
  ss<-dwt_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),80)
})
