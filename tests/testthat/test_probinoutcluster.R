library("testthat")
context("In/out cluster detection works correctly")


test_that("All parameters are correctly specified and trouble-shooted in clust_sim", {
    ##################
    #Set- up for test 
    ##################
    #set initial params
    n <- 208
    Time <- 5
    risk.ratio <- 1
    rMax <- 20
    nsim <- 2
    #create lassoresult
    thresh <- NULL
    #create rr
    rr = matrix(1, nrow=n, ncol=Time)
    rr[5:10,2:4] <- risk.ratio #create fake cluster
    
    #create fake coordinates
    set.seed(1)
    x <- rnorm(n, mean = 358.1, sd=33.8)
    y <- rnorm(n, mean = 4016, sd = 32.7)
    clusters <- clusters2df(x,y,rMax, utm=TRUE, length(x))
    n <- length(unique(clusters$center))
    numCenters <- n
    sparseMAT <- spacetimeMat(clusters, numCenters, Time)
    potcl <- dim(sparseMAT)[2]
    #create lassoresult objects
    select.qbic <- list(450, 500)
    select.qaic <- list(600, 625)
    select.qaicc <- list(600, 625)

    #create beta matrix #dims = c(66870, 450) and c(66870, 400)
    sparseMatrix(i = c(rep(51893,154),rep(51894,40), rep(66870,446)),
                         j = seq(1:640),
                         x = rnorm(640, mean = 0.03, sd = 0.04))
    a <- list(beta = sparseMatrix(i = c(rep(51893,154),rep(51894,40), rep(potcl,446)),
                                           j = seq(1:640),
                                           x = rnorm(640, mean = 0.03, sd = 0.04)))
    b <- list(beta = sparseMatrix(i = c(rep(51893,154),rep(51894,40), rep(potcl,446)),
                                  j = seq(1:640),
                                  x = rnorm(640, mean = 0.03, sd = 0.04)))
    lasso <- list(a,b)
    
    #put together into lassoresult
    lassoresult <- list(select.qbic = select.qbic, select.qaic=select.qaic, select.qaicc = select.qaicc,
        nsim = nsim, 
        lasso = lasso
    )
    
    #test the function - should not find anything inside and finds everything outside
    moo <- prob_inoutcluster(lassoresult, rr, risk.ratio, x, y, rMax, nsim,Time,thresh)
    expect_equal(moo$notinperc.bic, as.character(paste0(100,"%")))
    expect_equal(moo$notinperc.aic, as.character(paste0(100,"%")))
    expect_equal(moo$inperc.bic, as.character(paste0(0,"%")))
    expect_equal(moo$inperc.aic, as.character(paste0(0,"%")))
}
)
