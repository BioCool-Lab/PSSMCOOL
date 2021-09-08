context("checking that the length of pssm400 feature vector is equal to 400")

test_that("whether the pssm400 gives us the expected output",{
  ss<-pssm400(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),400)
})
