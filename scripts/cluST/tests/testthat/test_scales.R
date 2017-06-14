library("testthat")
context("Standardizing observed and expected")


test_that("Test simulation scales", {
    #set up
    Time <- 2
    theta = 1000
    nsim = 2
    period <- rep(seq(1,2),5)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(E, theta=1000)
    #set init
    init <- setVectors(period, expected, observed, Time, byrow=TRUE)
    ysim <- lapply(1:nsim, function(i) rnegbin(expected, theta = theta))
    
    #expectations
    expect_that(length(scale.sim(ysim, init, nsim, Time)), equals(nsim))
    
})


