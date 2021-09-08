context("checking that the length of SOMA_PSSM feature vector is equal to 160")

test_that("whether the SOMA_PSSM gives us the expected output",{
  ss<-SOMA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),160)
})
