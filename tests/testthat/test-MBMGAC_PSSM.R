context("checking that the length of MBMGACPSSM feature vectorS is equal to 560")

test_that("whether the MBMGACPSSM gives us the expected output",{
  ss<-MBMGACPSSM(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  expect_equal(length(ss),560)
})
