context("checking that the length of AC feature vector is equal to 340")

test_that("whether the AC function gives us the expected output",{
  ss<-pssm_ac(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),17)
  expect_equal(length(ss),340)
})
