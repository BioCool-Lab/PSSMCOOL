context("checking that the length of k_seperated_bigrame feature vector is equal to 400")

test_that("whether the k_seperated_bigrame function gives us the expected output",{
  ss<-k_seperated_bigrame(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),5)
  expect_equal(length(ss),400)
})
