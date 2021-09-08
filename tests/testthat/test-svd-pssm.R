context("checking that the length of SVD_PSSM feature vector is equal to 20")

test_that("whether the SVD_PSSM gives us the expected output",{
  ss<-SVD_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),20)
})
