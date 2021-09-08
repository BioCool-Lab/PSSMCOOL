context("checking that the length of CC feature vector is equal to 380")

test_that("whether the CC function gives us the expected output",{
  ss<-pssm_cc(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),19)
  expect_equal(length(ss),380)
})
