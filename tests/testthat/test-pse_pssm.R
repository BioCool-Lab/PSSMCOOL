context("checking that the length of pse_pssm feature vector is equal to 320")

test_that("whether the pse_pssm function gives us the expected output",{
  ss<-pse_pssm(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  expect_equal(length(ss),320)

})
