context("checking that the length of DPC_PSSM feature vector is equal to 420")

test_that("whether the DPC_PSSM function gives us the expected output",{
  ss<-DPC_PSSM(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  expect_equal(length(ss),420)
})
