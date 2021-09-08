context("checking that the length of trigrame_pssm feature vector is equal to 8000")

test_that("whether the trigrame_pssm gives us the expected output",{
  ss<-trigrame_pssm(system.file("extdata", "C7GSI6.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),8000)
})
