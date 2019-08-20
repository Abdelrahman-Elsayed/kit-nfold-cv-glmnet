

# install.packages("glmnet")
require(glmnet)
#####################################################################################
# Perform K iterations of n-fold cross-validation of a glmnet model fitting procedure

k.nfold.cv.glmnet=function(dset,           # a data.frame
                           yvar,           # name or numeric index of response variable in dset
                           xvars,          # names or numeric indices of predictor variables in dset
                           fam,            # indicates model family to fit, see help(glmnet) for more details
                           alpha=1,        # elastic net mixing parameter, alpha = 1 is lasso, alpha=0 is ridge, see help(glmnet) for more info
                           nfolds=10,      # number of cross-validation folds, default is 10, see help(cv.glmnet) for more details
                           k=1000)         # number of iterations of nfold cross validation, default is 1000
{
  # Check input data
  if (!is.data.frame(dset))
    stop("dset must be a data.frame")
  if (is.character(yvar))
  {
    if (!is.element(yvar,colnames(dset)))
      stop("yvar must be the numeric index or character name of a column of the response variable in dset.")
  }
  if (is.character(xvars))
  {
    if (any(!is.element(xvars,colnames(dset))))
      stop("xvars must be a vector of numeric indices or character names of the candidate predictor variables in dset.")
  }
  
  
  y=dset[,yvar]                   # extract response variable
  X=dset[,xvars]                  # extract candidate predictor variables
  X=as.matrix(X)                  # represent candidate predictors as a matrix
  B=matrix(NA,ncol(X),k)          # initialize matrix of coefficient estimates
  rownames(B)=colnames(X)         # label rows of coefficient estimates with variable names
  colnames(B)=paste0("CV_",1:k)   # label columns of coefficient estimates with 
  for (i in 1:k)                  # loop over k iterations of nfold cross-validation
  {
    message(paste0("Performing iteration ",i," of ",k, 
                   " iterations of ",nfolds,
                   "-fold cross-validation: ",
                   date()))
    fit.res=cv.glmnet(y=y,x=X,                  # perform the n-fold cross-validation of glmnet fit
                      nfolds=nfolds,            
                      family=fam)
    fit.est=as.numeric(coef(fit.res,            # extract the estimates from this iteration
                            s="lambda.min"))
    B[,i]=fit.est                               # place those estimates in the coefficient estimate matrix
  }
  
  class(B)="k.nfold.cv.glmnet"                  # assign a class to the result matrix for printing
  
  return(B)                                     # return the result matrix
}  # end of function k.nfold.cv.glmnet


#################################################################
# Display results of k iterations of n-fold cross-validation
print.k.nfold.cv.glmnet=function(B,alpha=0.05,
                                 max.rows=10)# B is the result of k.nfold.cv.glmnet

{
  non.zero=rowMeans(B!=0)
  mean.coef=rowMeans(B)
  median.coef=apply(B,1,median)
  q1.coef=apply(B,1,quantile,0.25)
  q3.coef=apply(B,1,quantile,0.75)
  lower.bound=apply(B,1,quantile,alpha/2)
  upper.bound=apply(B,1,quantile,1-alpha/2)
  
  
  res=cbind(non.zero=non.zero,
            mean.coef=mean.coef,
            median.coef=median.coef,
            q1.coef=q1.coef,
            q3.coef=q3.coef,
            lower.bound=lower.bound,
            upper.bound=upper.bound)
  colnames(res)[6:7]=paste0(c("lower.","upper."),
                            1-alpha,".bound")
            
  ord=rev(order(res[,"non.zero"],abs(res[,"mean.coef"])))
  res=res[ord,]
  n.non.zero=sum(non.zero!=0)
  if (n.non.zero>max.rows)
    warning(paste0("There are ",n.non.zero,
                   " variables with a non-zero ",
                   " coefficient in at least one iteration.  ",
                   "Only the top ",max.rows," are shown here."))
  print(head(res,max.rows))
}

##########################################
# Simulated Example
n=200                              # number of patients
g=50                               # number of genes
X=matrix(rnorm(n*g),n,g)           # expression matrix
colnames(X)=paste0("gene_",1:g)    # label the columns
true.stime=rexp(n,exp(2*X[,1]))    # simulated true survival times depend on gene 1
fu.time=runif(n,3,8)               # simulated follow-up times
obs.stime=pmin(fu.time,true.stime) # observed survival time is minimum of follow-up time and true survival time
obs.status=as.numeric(true.stime<=fu.time)  # observed status

library(survival)
obs.surv=Surv(obs.stime,obs.status)  
dset=cbind.data.frame(obs.surv=obs.surv,X)

sim.res=k.nfold.cv.glmnet(dset,
                          yvar="obs.surv",
                          xvars=grep("gene",colnames(dset)),
                          fam="cox",k=100)
sim.res
