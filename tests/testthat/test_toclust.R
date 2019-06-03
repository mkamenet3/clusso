library("testthat")
library("MASS")
context("Testing toclusso")



test_that("Length of inputs the same",{
    set.seed(2)
    period <- rep(seq(1,2),5)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    Time = 2
    dat <- cbind.data.frame(expected, observed, period)
    
    expect_warning(toclusso(dat,expected, observed, period))
    expect_error(toclusso(dat,expected, observed[-1], as.factor(period)))
    expect_error(toclusso(c(expected,9), observed, as.factor(period)))
    expect_error(toclusso(C(expected,9), observed[-1], as.factor(period)))
})



