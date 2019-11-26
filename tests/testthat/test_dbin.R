library("testthat")
context("Binomial log-likelihood")

test_that("Binomial log-lik and phat",{
    #setup
    set.seed(2)
    period <- c(rep("1",5),rep("2",5))
    observed <- MASS::rnegbin(n = length(period), mu=20, theta=2)
    pop <-  MASS::rnegbin(n = length(period), mu=200, theta=2)
    ids <- rep(1:5, times=2)
    m <- glm(cbind(observed, (pop-observed)) ~ period, family="binomial")
    phats <- pihat(model.matrix((m))%*%coef(m))
    dbin(observed, pop,phats)
    expect_equal(dbin(observed, pop,phats), -68.6853, tolerance = 0.001)
})
