context("checking that the dimention of dataframe obtained from disulfid function is 45*521 ")

test_that("whether the disulfid gives us the expected output",{
  ss<-disulfid(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(dim(ss),c(28,521))
})
