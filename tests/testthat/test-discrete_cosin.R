context("checking that the length of Discrete_Cosine_Transform feature vector is equal to 400")

test_that("whether the Discrete_Cosine_Transform gives us the expected output",{
  ss<-Discrete_Cosine_Transform(system.file("extdata", "C7GQS7.txt.pssm", package="PSSMCOOL"))
  expect_equal(length(ss),400)
})
