context("checking that the length of PSSM_SD and PSSM_SEG feature vectors are equal to 80,100 respectivly")

test_that("whether the PSSM_SD,PSSM_SEG functions gives us the expected output",{
  ss<-PSSM_SD(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"))
  ss2<-pssm_seg(system.file("extdata","C7GQS7.txt.pssm",package="PSSMCOOL"),3)
  expect_equal(length(ss[[2]]),80)
  expect_equal(length(ss2),100)
})
