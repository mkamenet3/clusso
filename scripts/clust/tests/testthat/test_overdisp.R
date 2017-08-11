library("testthat")
context("Estimating overdispersion for Quasi-Poisson model")

test_that("Test overdispersion model is from glm",{
    #setup
    period <- c(rep("1",5),rep("2",5))
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    m <- lm(observed ~ 1 + as.factor(period))
    expect_error(overdisp(m))
})
