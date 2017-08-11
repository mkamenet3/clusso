library("testthat")
context("Standardizing observed and expected")


test_that("Test simulation scales", {
    #set up
    set.seed(2)
    Time <- 2
    theta = 1000
    nsim = 2
    period <- rep(seq(1,2),5)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    #set init
    init <- setVectors(period, expected, observed, Time, byrow=TRUE)
    ysim <- lapply(1:nsim, function(i) rnegbin(expected, theta = theta))
    ysim2 <- ysim
    ysim2[[1]] <- ysim[[1]][-1]
    #expectations
    expect_that(length(scale.sim(ysim, init, nsim, Time)), equals(nsim))
    expect_error(scale.sim(ysim2, init, nsim, Time))
})

test_that("Test scales (non-sim)", {
    #set up
    set.seed(2)
    Time <- 2
    theta = 1000
    nsim = 2
    period <- rep(seq(1,2),5)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    #set init
    init <- setVectors(period, expected, observed, Time, byrow=TRUE)
    #expectations
    expect_error(length(scale(init, 4)))
})


