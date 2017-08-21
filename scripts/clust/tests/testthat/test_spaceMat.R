library("testthat")
context("Creating sparse spatial matrix")

test_that("Input is actually a matrix", {
    #set-up
    rMax = 20
    #utm example
    x_utm <- c(399786.6, 360917.0, 385175.1, 371603.4, 388154.2, 375023.3)
    y_utm <- c(4047756, 4023885, 4025749, 4018172, 4047900, 4068053)
    #lat/long example
    x_latlon <- c(36.569996, 36.350001, 36.370002, 36.299997, 36.570002, 36.749997)
    y_latlon<- c(7.88, 7.45, 7.72, 7.57, 7.57, 7.60)
    clusterdf <- clusters2df(x_utm, y_utm, rMax, utm=TRUE, length(x_utm))
    is.sparseMatrix <- function(x) is(x, 'sparseMatrix') 
    
    #great expectations
    expect_true(is.sparseMatrix(spaceMat(clusterdf, length(x_utm))))
})