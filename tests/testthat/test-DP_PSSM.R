context("checking that the length of DP_PSSM feature vector is equal to 120")

test_that("whether the DP_PSSM function gives us the expected output",{
  ss<-DP_PSSM(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  expect_equal(length(ss),120)
})
