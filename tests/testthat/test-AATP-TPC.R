context("checking that the length of AATP_TPCC feature vector is equal to 420")

test_that("whether the AATP_TPCC gives us the expected output",{
  ss<-AATP_TPCC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),420)
})
