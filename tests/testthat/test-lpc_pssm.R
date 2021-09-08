context("checking that the length of LPC_PSSM feature vector is equal to 280")

test_that("whether the LPC_PSSM gives us the expected output",{
  ss<-LPC_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),280)
})
