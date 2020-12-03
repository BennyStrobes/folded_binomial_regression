require(rstan)
# rstan_options(auto_write = TRUE)


FOLDED_BINOMIAL_GLM=stan_model(file="folded_binomial_glm.stan", save_dso=F, auto_write=F)


# logistic (sigmoid) function
logistic=function(g) 1/(1+exp(-g))

# ys: numerator counts
# ns: denominator counts
# xFull: matrix of covariates for full model. First column must be ones. 
# xNull: matrix of covariates for null model. First column must be ones. 
# Prior on concentration parameter is Gamma(concShape,concRate)
foldedBinomialGLM=function(ys,ns,xFull,xNull) {
  stopifnot(all(xNull==xFull[,1:ncol(xNull)]))
  stopifnot(all(xNull[,1]==1))

  
  dat=list(N=length(ys),P=ncol(xNull),ys=ys,ns=ns,x=xNull)
  
  # Fit null model
  fit_null <- optimizing(FOLDED_BINOMIAL_GLM, data=dat, algorithm="BFGS", hessian=T, as_vector=F)
  
  # Initialize alternative model using null model
  betaInit=numeric(ncol(xFull))
  betaInit[1:ncol(xNull)]=fit_null$par$beta
  initFull=list(beta=betaInit)
  
  # Fit alternative model
  datFull=dat
  datFull$x=xFull
  datFull$P=ncol(xFull)
  fit_full <- optimizing(FOLDED_BINOMIAL_GLM, data=datFull, init=initFull, algorithm="BFGS", hessian=T, as_vector=F)
  
  loglr=fit_full$value - fit_null$value
  
  list( loglr=loglr, lrtp=pchisq( 2.0*loglr, lower.tail = F , df=ncol(xFull)-ncol(xNull) ), fit_full=fit_full, fit_null=fit_null )
}

# Function to extract coefficients, standard errors, and Wald p-values
get_glm_coefs=function(res) {
  p=length(res$par$beta)
  variance=robustSolve(-res$hessian)
  dimnames(variance)=dimnames(res$hessian)
  betase=sqrt(diag(variance))[paste0("beta.",1:p)]
  beta=res$par[paste0("beta[",1:p,"]")]
  zscore=res$par$beta / betase
  data.frame(co=res$par$beta,se=betase,p=2.0*pnorm(-abs(zscore)))
}

null_test_foldedbinomial_glm=function() {
  nsamp=100
  ns=rpois(nsamp,lambda=10)
  ys=rbinom(nsamp,ns,logistic(.3+.5*rnorm(nsamp)))
  ys_folded = pmin(ys, ns-ys)
  xNull=cbind(numeric(nsamp)+1)
  xFull=cbind(xNull,runif(nsamp))
  anova_result=foldedBinomialGLM(ys_folded,ns,xFull,xNull)
  print(anova_result$lrtp)
  #get_glm_coefs(anova_result$fit_null)
  #get_glm_coefs(anova_result$fit_full)
}
test_foldedbinomial_glm=function() {
  nsamp=100
  ns=rpois(nsamp,lambda=10)
  env = rnorm(nsamp)
  ys=rbinom(nsamp,ns,logistic(.3+.5*env))
  ys_folded = pmin(ys, ns-ys)
  xNull=cbind(numeric(nsamp)+1)
  xFull=cbind(xNull,env)
  anova_result=foldedBinomialGLM(ys_folded,ns,xFull,xNull)
  print(anova_result$lrtp)
  #get_glm_coefs(anova_result$fit_null)
  #get_glm_coefs(anova_result$fit_full)
}

null_test_foldedbinomial_glm()
test_foldedbinomial_glm()
