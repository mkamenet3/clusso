library("testthat")
library("MASS")
context("Poisson log-likelihood")


test_that("inputs are non-zero",{
    #expect_that(round(dpoisson(c(5,10),c(1,1),c(1,5)),3), equals(10.094))
    expect_error(round(dpoisson(c(5,10),c(0,1),c(1,5)),3))
    #expect_that(round(dpoisson(5,1,5),3), equals(3.047))
    expect_warning(dpoisson(5,1,0))
})



