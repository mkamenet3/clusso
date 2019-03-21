

clustplotlasso <- function(){
    
    #plot(resreal$lassoresult.p.s$lasso_out$lambdas, resreal$lassoresult.p.s$lasso_out$coefs_bic)
    loglambdas <- log(resreal$lassoresult.p.s$lasso_out$lambdas)
    
    plot(log(resreal$lassoresult.p.s$lasso_out$lambdas),resreal$lassoresult.p.s$lasso_out$coefs_bic[,1], type="l")
    lines(add=TRUE)
    
    plot(resreal$lassoresult.p.st$lasso,xvar = "lambda",  main="test");abline(v = loglambdas,lty=1, lwd=3)
    
    plot_glmnet(resreal$lassoresult.p.st$lasso)
    
    plot_glmnet(resreal$lassoresult.p.st$lasso, xvar="rlambda");abline(v = loglambdas,lty=1, lwd=3)
}