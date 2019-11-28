library("testthat")
library("MASS")
context("Setting observed, expected, and time period vectors")



test_that("Data imported correct way (byrow check)",{
    set.seed(2)
    period1 <- c(rep("1",5),rep("2",5))
    period2 <- rep(seq(1,2),5)
    expected <- MASS::rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- MASS::rnegbin(expected, theta=1000)
    Time = 2
    covars <- NULL
    expect_error(setVectors(period1, expected, observed, covars, Time, byrow=TRUE))
})




