library("testthat")
library("MASS")
context("Testing toclust")



test_that("Length of inputs the same",{
    set.seed(2)
    period <- rep(seq(1,2),5)[-1]
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    Time = 2
    expect_error(toclust(expected, observed, period))
    expect_error(toclust(expected, observed[-1], period))
    expect_error(toclust(c(expected,9), observed, period))
    expect_error(toclust(C(expected,9), observed[-1], period))
})



