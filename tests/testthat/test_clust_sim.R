library("testthat")
context("Simulation of cluster lasso")

test_that("All parameters are correctly specified and trouble-shooted in clust_sim", {
    #set up
    x <- c(399786.6, 360917.0, 385175.1, 371603.4, 388154.2, 375023.3,383905.5, 392283.1, 412155.6, 412654.0)
    y <- c(4047756, 4023885, 4025749, 4018172, 4047900, 4068053, 4064599, 4020110, 4032090, 4080899)
    Time <- 2
    theta = 1000
    nsim = 2
    center = 1
    radius = 10
    overdispfloor = TRUE
    risk.ratio = 2
    threshold = 0.9
    period <- rep(seq(1,2),5)
    timeperiod <- 1:2
    covars <- NULL
    r.Max <- 12
    theta <- 1000
    #expected <- lapply(1:nsim, function(i) rnegbin(n = 10,mu = 15,theta = 1000))
    #observed <- lapply(1:nsim, function(i) rnegbin(n = 10,mu = 20,theta = 1000))
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(n = 10,mu = 20,theta = 1000)
    YSIM <- lapply(1:nsim, function(i) rnegbin(n = 10,mu = 15,theta = 1000))
    init <- setVectors(period, expected, observed, covars,  Time, byrow=TRUE)
    df <- cbind.data.frame(expected = expected, observed = observed, timeperiods = period)
    clst <- toclust(df, expected = df$expected, observed = df$observed, timeperiod = df$timeperiods, covars = FALSE)
    
    expect_error(clust_sim(x_utm, y_utm, r.Max, Time, nsim, center, radius, risk.ratio, timeperiod, utm=TRUE, byrow=TRUE,
                           threshold = 0.5,
                           space = "both", theta, nullmod = NULL, overdispfloor=FALSE))
})


test_that("Vectors_space_sim correctly collapses from space-time to space.", {
    #set up
    set.seed(2)
    Time <- 2
    theta = 1000
    nsim = 2
    period <- rep(seq(1,2),5)
    expected <- list(rnegbin(n = 10,mu = 15,theta = 1000))
    observed <- rnegbin(unlist(expected), theta=1000)
    x <- c(399786.6, 360917.0, 385175.1, 371603.4, 388154.2)
    Ex <- list(unlist(expected), unlist(expected))
    covars.s <- NULL
    #set init
    init <- setVectors(period, unlist(expected), observed, covars.s, Time, byrow=TRUE)
    ysim <- lapply(1:nsim, function(i) rnegbin(unlist(expected), theta = theta))
    
    expect_equal(length(vectors_space_sim(x,Ex,ysim, Time,init)),2)
})
