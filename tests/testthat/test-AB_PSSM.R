context("checking that the length of AB-PSSM feature vectorS is equal to 400")

test_that("whether the AB-PSSM gives us the expected output",{
  ss<-AB_PSSM(system.file("extdata","C7GRQ3.txt.pssm",package="PSSMCOOL"))
  expect_equal(length(ss[[1]]),400)
  expect_equal(length(ss[[2]]),400)
})
