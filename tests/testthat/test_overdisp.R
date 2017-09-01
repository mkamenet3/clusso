library("testthat")
context("Estimating overdispersion for Quasi-Poisson model")

test_that("Test overdispersion model is from glm",{
    #setup
    set.seed(2)
    period <- c(rep("1",5),rep("2",5))
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    m <- glm(observed ~ 1 + as.factor(period),family = "poisson")
    expect_error(overdisp(m))
})

test_that("Test underdispersion is handled properly with overdispfloor argument",{
    #setup
    set.seed(2)
    period <- c(rep("1",5),rep("2",5))
    expected <- rnegbin(n = 10,mu = 0.05,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    m <- glm(observed ~ 1 + as.factor(period), family="poisson")
    expect_message(overdisp(m, sim=FALSE))
    expect_message(overdisp(m, sim=FALSE, overdispfloor=FALSE))
})


test_that("Test sim work",{
    #setup
    set.seed(2)
    period <- c(rep("1",5),rep("2",5))
    expected <- list(rnegbin(n = 10,mu = 15,theta = 1000),
                     rnegbin(n = 10,mu = 15,theta = 1000))
    observed <- list(rnegbin(expected[[1]], theta=1000),
                     rnegbin(expected[[1]], theta=1000))
    m <- lapply(1:nsim, function(i) glm(observed[[i]] ~ 1 + as.factor(period)))
    expect_gt(overdisp(m), 1)
    expect_gt(overdisp(m, sim=TRUE), 1)
})


test_that("Test sim works with overdispfloor",{
    #setup
    set.seed(2)
    period <- c(rep("1",5),rep("2",5))
    expected <- list(rnegbin(n = 10,mu = 1,theta = 1000),
                     rnegbin(n = 10,mu = 1,theta = 1000))
    observed <- list(rnegbin(expected[[1]], theta=1000),
                     rnegbin(expected[[2]], theta=1000))
    m <- lapply(1:nsim, function(i) glm(observed[[i]] ~ 1 + as.factor(period)))
    expect_message(overdisp(m, sim=TRUE, overdispfloor=TRUE))
    expect_message(overdisp(m, sim=TRUE))
})