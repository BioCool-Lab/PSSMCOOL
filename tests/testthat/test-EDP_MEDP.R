context("checking that the length of EDP_MEDP feature vectorS is equal to 20,400,420,respectively")

test_that("whether the EDP_MEDP gives us the expected output",{
  ss<-EDP_MEDP(paste0(system.file("extdata",package="PSSMCOOL"),"/C7GS61.txt.pssm"))
  expect_equal(length(ss),420)
})

