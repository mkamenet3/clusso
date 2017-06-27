library("testthat")
context("Simulation of cluster lasso")

test_that("All parameters are correctly specified and trouble-shooted in clust.sim.all", {
    #set up
    x <- c(399786.6, 360917.0, 385175.1, 371603.4, 388154.2, 375023.3,383905.5, 392283.1, 412155.6, 412654.0)
    y <- c(4047756, 4023885, 4025749, 4018172, 4047900, 4068053, 4064599, 4020110, 4032090, 4080899)
    Time <- 2
    theta = 1000
    nsim = 2
    center = 1
    radius = 10
    risk.ratio = 2
    threshold = 0.9
    period <- rep(seq(1,2),5)
    timeperiod <- 1:2
    expected <- lapply(1:nsim, function(i) rnegbin(n = 10,mu = 15,theta = 1000))
    observed <- lapply(1:nsim, function(i) rnegbin(n = 10,mu = 20,theta = 1000))
    YSIM <- lapply(1:nsim, function(i) rnegbin(n = 10,mu = 15,theta = 1000))
    
    
    expect_error(clust.sim.all(x_utm, y_utm, period, expected, observed, Time, nsim, center,
                               radius, risk.ratio, timeperiod, utm=TRUE, byrow=TRUE, threshold))
    expect_error(vectors.space(x,expected,Time, timeperiod))
    
})
