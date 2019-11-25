library("testthat")
context("Calculating expected counts")

test_that("Dataframe with expected counts",{
    #setup
    set.seed(2)
    period <- c(rep("1",5),rep("2",5))
    observed <- MASS::rnegbin(n = length(period), mu=20, theta=2)
    pop <-  MASS::rnegbin(n = length(period), mu=200, theta=2)
    ids <- rep(1:5, times=2)
    expect_is(calcexpected(observed, pop, period, ids), "data.frame")
})