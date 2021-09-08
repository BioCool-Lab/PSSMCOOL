context("checking that the length of single_Average feature vector is equal to 400")

test_that("whether the single_Average gives us the expected output",{
  ss<-single_Average(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),400)
})
