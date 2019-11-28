library("testthat")
context("Standardizing observed and expected")

test_that("Test scales", {
    #set up
    set.seed(2)
    Time <- 2
    theta = 1000
    nsim = 2
    period <- rep(seq(1,2),5)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    covars <- NULL
    #set init
    init <- setVectors(period, expected, observed, covars, Time, byrow=TRUE)
    #expectations
    expect_error(length(scale(init, 4)))
})


