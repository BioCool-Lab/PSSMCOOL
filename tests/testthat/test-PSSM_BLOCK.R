context("checking that the length of PSSMBLOCK feature vector is equal to 300")

test_that("whether the PSSMBLOCK function gives us the expected output",{
  ss<-PSSMBLOCK(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),5)
  ss2<-Averag_Block(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  expect_equal(length(ss),300)
  expect_equal(length(ss2),400)
})
