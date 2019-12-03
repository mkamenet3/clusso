library("testthat")
library("MASS")
context("Testing toclusso")



test_that("Length of inputs the same",{
    set.seed(2)
    period <- rep(seq(1,2),5)
    periodF <- as.factor(period)
    expected <- rnegbin(n = 10,mu = 15,theta = 1000)
    observed <- rnegbin(expected, theta=1000)
    covars <- FALSE
    id <- NULL
    requiredcolNames <- c("expected", "observed", "period")
    requiredcolNamesF <- c("expected", "observed", "periodF")
    Time = 2
    dat <- cbind.data.frame(expected, observed, period, periodF)
    
    expect_warning(toclusso(dat,expected, observed, period, covars, id, requiredcolNames))
    expect_error(toclusso(dat,expected, observed[-1], periodF, covars, id, requiredcolNamesF))
    expect_error(toclusso(c(expected,9), observed, periodF, covars, id, requiredcolNamesF))
    expect_error(toclusso(C(expected,9), observed[-1], periodF, covars, id, requiredcolNamesF))
})



