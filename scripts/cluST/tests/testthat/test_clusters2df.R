library("testthat")
library("MASS")
context("Creating dataframe of potential clusters based on predefined max radius")


test_that("Warning that coordinates are actually UTM or Lat/Long", {
    rMax = 20
    #utm example
    x_utm <- c(399786.6, 360917.0, 385175.1, 371603.4, 388154.2, 375023.3)
    y_utm <- c(4047756, 4023885, 4025749, 4018172, 4047900, 4068053)
    #lat/long example
    x_latlon <- c(36.569996, 36.350001, 36.370002, 36.299997, 36.570002, 36.749997)
    y_latlon<- c(7.88, 7.45, 7.72, 7.57, 7.57, 7.60)
    
    #test
    expect_error(clusters2df(x_utm, y_utm, rMax, utm=FALSE, length(x_utm)))
    expect_error(clusters2df(x_latlon, y_latlon, rMax, utm=TRUE, length(x_latlon)))
    
})




