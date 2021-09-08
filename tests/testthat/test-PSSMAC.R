context("checking that the length of PSSMAC feature vector is equal to 200")

test_that("whether the PSSMAC gives us the expected output",{
  w<-PSSMAC(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(w),200)
})
