context("checking that the length of rpssm feature vector is equal to 110")

test_that("whether the rpssm gives us the expected output",{
  ss<-rpssm(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),110)
})
