context("checking that the length of DMACA_PSSM feature vectorS is equal to 210")

test_that("whether the DFMCA_PSSM gives us the expected output",{
  ss<-DFMCA_PSSM(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"),4)
  expect_equal(length(ss),210)
})
