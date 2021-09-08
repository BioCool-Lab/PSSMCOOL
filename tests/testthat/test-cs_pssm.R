context("checking that the length of CS_PSe_PSSM feature vector is equal to 40")

test_that("whether the CS_PSe_PSSM function gives us the expected output",{
  ss<-CS_PSe_PSSM(system.file("extdata","C7GSI6.txt.pssm",package="PSSMCOOL"),"cspssm")
  expect_equal(length(ss),40)
})
