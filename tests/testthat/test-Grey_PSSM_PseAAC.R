context("checking that the length of grey_PSSM feature vector is equal to 100")

test_that("whether the grey_PSSM function gives us the expected output",{
  ss<-grey_pssm_pseAAC(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  expect_equal(length(ss),100)
})
