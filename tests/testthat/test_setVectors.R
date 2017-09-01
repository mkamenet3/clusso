library("testthat")
library("MASS")
context("Setting observed, expected, and time period vectors")



test_that("Data imported correct way (byrow check)",{
    period1 <- c(rep("1",5),rep("2",5))
    period2 <- rep(seq(1,2),5)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    Time = 2
    covars <- NULL
    expect_warning(setVectors(period1, expected, observed, covars, Time, byrow=TRUE))
    expect_warning(setVectors(period2, expected, observed, covars, Time, byrow=FALSE))
})




