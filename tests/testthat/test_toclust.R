library("testthat")
library("MASS")
context("Testing toclust")



test_that("Length of inputs the same",{
    set.seed(2)
    period <- rep(seq(1,2),5)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    Time = 2
    dat <- cbind.data.frame(expected, observed, period)
    
    expect_warning(toclust(dat,expected, observed, period))
    expect_error(toclust(dat,expected, observed[-1], as.factor(period)))
    expect_error(toclust(c(expected,9), observed, as.factor(period)))
    expect_error(toclust(C(expected,9), observed[-1], as.factor(period)))
})



