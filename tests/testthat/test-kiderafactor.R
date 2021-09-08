context("checking that the length of kiderafactor feature vector is equal to 80")

test_that("whether the kiderafactor gives us the expected output",{
  ss<-kiderafactor(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),c(2,3,8,9))
  expect_equal(dim(ss)[2],170)
})
