##############################################
#
# Functions for ensemble end-member mixing analysis:
#
#     eemma       performs ensemble end-member mixing analysis
#
#     eiv         implements a modified version of Fuller's procedure for correcting biases in multiple regression coefficients
#                     that arise due to errors in the explanatory variables ("x-axis" errors)
#
#     ols         an ordinary least squares multiple regression routine, which can be used instead of eiv when the errors/uncertainties 
#                     in the explanatory variables are a small fraction of their variances, and thus no "errors-in-variables" correction is needed




# version 1.0  build 2023.03.27
# Author: James Kirchner, ETH Zurich
#
# Copyright (C) 2023 ETH Zurich and James Kirchner
# Public use of this script is permitted under GNU General Public License 3 (GPL3); for details see <https://www.gnu.org/licenses/>
# 
#
# READ THIS CAREFULLY:
# ETH Zurich and James Kirchner make ABSOLUTELY NO WARRANTIES OF ANY KIND, including NO WARRANTIES, expressed or implied, that this software is
#    free of errors or is suitable for any particular purpose.  Users are solely responsible for determining the suitability and
#    reliability of this software for their own purposes.


########################################################################################
# ALSO READ THIS:
#
# This code implements Ensemble End-Member Mixing Analysis, as described in:
# Kirchner, J.W. (2023), Mixing models with multiple, overlapping, or incomplete end-members, 
# quantified using time series of a single tracer, Geophysical Research Letters. 
#
# Users should cite that paper for the concepts and mathematics.
#
# They should also cite Kirchner, J.W. (2023), eemma, an R script for ensemble end-member mixing analysis, www.envidat.ch/datset/eemma, for the code itself.
########################################################################################




library(dplyr)  # if dplyr is not already installed, do install.packages("dplyr")





################################################
# eiv: performs errors-in-variables regression #
################################################
eiv <- function(
    y,                                  # vector of response variable
    X,                                  # matrix of one or more explanatory variables
    no.intercept = FALSE,               # flag for whether regression should be without intercept
    err.all = 1,                        # measurement error in x's and y, expressed as standard deviation
    err.y = err.all,                    # measurement error in y, expressed as standard deviation
    err.x = rep(err.all, NCOL(X)),      # vector of measurement errors in x's, expressed as standard deviation
    rxx = diag(NCOL(X)),                # correlation matrix of measurement errors in x's (correlations, not covariances!)
    rxy = rep(0, NCOL(X)),              # correlation between measurement errors in y and each of the x's (correlations, not covariances!)
    Szz = NULL,                         # covariance matrix of Z=cbind(y,X) (for calls from eemma)
    uncert = c("relative", "absolute"), # are uncertainties absolute, or relative?
    alpha = 1                           # parameter with value between 1 and 2, inclusive
)

{ 
  # Performs "errors in variables" linear regression, correcting for biases introduced by errors in explanatory variables.
  # Implements Fuller's bias correction procedure using equations 2.5.13 and 2.5.14 from pages 173 and 174 of Fuller, Measurement Error Models (1987),
  # modified as detailed below.
  
  # receives as input:
  #
  #     y             vector containing the response variable (n rows, 1 column)
  #
  #     X             matrix, data table, or data frame containing one or more explanatory variables (n rows, k columns)
  #                       Column names, if available, are retained for output; otherwise names "x1", "x2", etc. are assigned.
  #
  #                   y and X must have the same number of rows.  NA values are allowed, but only complete cases will be analyzed.
  #
  #     no.intercept  if TRUE, no intercept is fitted.  Removing the intercept can help in constraining the coefficient of (maximum 1!)
  #                       explanatory variable that is constant or nearly constant.
  #
  #     err.all       standard deviation of measurement/sampling errors in x's and y.  If only err.all is provided, it is assumed to 
  #                       apply to both y and all the x's.  Default value is 1, with (by default) implicit rescaling as invoked by 
  #                       the default uncert=="relative" flag, below.
  #
  #     err.y         standard deviation of measurement/sampling error in y.  If provided, supersedes err.all.
  #
  #     err.x         standard deviation(s) of measurement/sampling errors in x's.  If provided, supersedes err.all.
  #                       Must be a vector of length NCOL(X), or a single value.
  #                       If a single value rather than a vector is supplied, it will be copied for all of the x's.
  #
  #     rxx           correlation matrix of measurement/sampling errors in x's (k by k matrix). 
  #                       If rxx is not specified, all cross-correlations are assumed to be zero.
  #
  #     rxy           vector of correlations between measurement/sampling errors between y and each of the x's.
  #                       If rxy is not specified, these cross-correlations are assumed to be zero.
  # 
  #     Szz           optional error covariance matrix for Z=cbind(y,X).  If provided, supersedes err.x, err.y, rxx, and rxy.
  #                       Intended primarily to simplify calling from eemma.
  #
  #     uncert        flag for whether the relative, or absolute, values of the uncertainties have been specified.  If uncert=="relative", 
  #                       uncertainties are expressed as relative values (for example, if uncertainty in the mixture measurements is twice the 
  #                       uncertainty in the end-member measurements, err.x and err.y can be any numbers in the ratio of 1:2) and these are
  #                       then implicitly rescaled in the fitting procedure.
  #                       Can be abbreviated.  Default is "relative".  
  #
  #     alpha         alpha parameter for Fuller's method.  Must be non-negative. 
  #                       Lower values give smaller bias but larger parameter variance; higher values give more bias but smaller parameter variance.
  #                       Values between 1 and 2, inclusive, have performed well in benchmark tests.
  
  
  # returns a list as output, with the following elements:
  #
  #     beta          EIV regression coefficients
  #
  #     se            standard errors of EIV regression coefficients
  #
  #     s2e           residual variance of EIV regression 
  #
  #     beta_cov      covariance matrix of regression coefficients
  #
  #     e             vector of residuals
  #
  #     lambda        lambda coefficient for Fuller's method
  
  
  
  
  #########################
  # root-finding function #
  #########################
  find_root <- function(cov.mat, err.mat) {
    # finds approximate smallest root lambda of det(cov.mat - lambda*err.mat), if there is one, between 0 and ulim, 
    # where:
    # cov.mat is a covariance matrix of measurements
    # err.mat is a the covariance matrix of the measurement/sampling errors in these vectors
    # ulim is the upper limit of the search domain (but if a solution is not found there, we search higher)
    
    func <- function(x) det(cov.mat - x*err.mat)
    
    npts <- 100 # number of points we will search for
    
    # the upper limit of our possible lambda values will be the one that would create
    # non-positive values anywhere on the diagonal of the matrix cov.mat-lambda*err.mat, because these are variances
    
    ulim <- 1.0/max(diag(err.mat)/diag(cov.mat), na.rm=TRUE)  # this is the upper limit of lambda 
    # I'm doing it this way b/c some elements of diag(err.mat) may be zero
    
    if (ulim==0) return(1) # if all the error variances are 0, return lambda=1 (since it will have no effect anyway)
    
    # here we adopt a brute-force approach to searching for roots
    x <- seq(from=0, to=ulim*1.02, length.out=npts)  # start with a mesh of npts in the interval of interest
    f <- sapply(x, func)  # calculate function values at each of these points
    
    # now search for zero-crossing roots of func
    if (min(f)<0) {
      i <- min(which(f<0))        # find the index of the first value below zero
      root <- uniroot(func, lower=x[i-1], upper=x[i])$root  # and find the root by Brent's method
    } else {
      root <- ulim                # if there is no such root, assign the upper limit (since any root must be less than this)
    }
    
    # now search for asymptotic roots of func (where it may parabolically approach zero but does not cross it)
    # we define an asymptotic root as any point that is less than 1% of the x=0 value (which should be the largest f for x values below the first root)
    fthresh <- max(0.01*f[1], min(f[f>0]))                                # or if the minimum positive value of f is bigger, use that as the threshold
    localmin <- ((f<=fthresh)&(f<dplyr::lag(f))&(f<dplyr::lead(f)))       # find the index of all x's that are local minima and less than this threshold
    if (sum(localmin, na.rm=TRUE)>0) localmin <- x[min(which(localmin))]  # if there are any local minima, pick the first one
    else localmin <- ulim
    
    if (localmin < root) return(localmin) else return(root)   # return the smallest of the zero-crossing roots or the asymptotic roots
    
  }  
  #################    
  # end find_root #
  #################
  
  
  
  
  
  ##############
  # check inputs
  
  uncert <- match.arg(uncert)  # uncertainties relative or absolute?
  
  y <- as.matrix(y)            # make these matrices in case they were provided as vectors, data frames, etc.
  X <- as.matrix(X)
  
  if (NROW(y)!=NROW(X)) stop("fatal error in eiv: X and y have different n")
  
  k <- NCOL(X)                  # number of X variables in our problem
  
  if(is.null(colnames(X))) xnames <- paste0("x",1:k)   # copy names of x variables if they are given; otherwise name them x1, x2, etc.
  else xnames <- colnames(X)
  
  # now cbind X and y to make Z, and retain only complete cases
  Z <- cbind(y, X)
  y <- y[complete.cases(Z),]
  X <- X[complete.cases(Z),]
  Z <- Z[complete.cases(Z),]
  
  n <- sum(complete.cases(Z))
  if (n<(k+1)) stop("Fatal error in eiv: too few complete cases -- no degrees of freedom left!")
  
  if (length(err.x)==1) err.x <- rep(err.x, NCOL(X))    # if err.x was specified as a single value, copy it for all the X variables
  else if (length(err.x)!=NCOL(X)) stop("Fatal error in eiv: err.x must be a single value, or a vector of the same length as NCOL(X)")
  
  
  #############################################
  # add dummy variable for intercept if desired
  
  if (no.intercept==FALSE) {
    k <- k+1                                 # the width of our design matrix will now expand by 1
    X <- as.matrix(cbind(X,1))               # add a column of 1's to the design matrix X
    rxy <- rbind(rxy, 0)                     # add a zero to Sxy
    rxx <- rbind(cbind(rxx,0),0)             # add a row and column of zeroes to Sxx
    rxx[k,k] <- 1                            # except we need a 1 at the lower right cell of Sxx
    err.x <- c(err.x, 0)                     # add a zero to the uncertainty vector of the x's
    if (!is.null(Szz)) Szz <- rbind(cbind(Szz,0),0)  # add a row and column of zeroes to Szz, if it is specified
    xnames <- c(xnames, "intercept")         # and lastly, add "intercept" to the column names
    Z <- cbind(y, X)                         # need to rebuild Z
  }
  
  
  ###############################
  # calculate covariances and mean crossproducts
  
  # Note that mZZ is a *covariance* and not a mean square, as it is shown in sect. 2.5.2 of Fuller.
  # A mean square will not work here, if the means of the x's are not all zero.
  mZZ <- cov(Z)*(n-1)/n   # rescaled for consistency with MXX and MXY
  MXX <- crossprod(X)/n
  MXY <- crossprod(X,y)/n
  
  
  ####################################################
  # now compile error covariance matrix Szz (covariance matrix of Z; Fuller calls this Saa, See, or gamma_ee)
  
  if (is.null(Szz)) {
    if (det(rxx)<=0) stop("Fatal error: error correlation matrix must be positive definite!")
    if (any(diag(rxx)!=1)) stop("Diagonal elements of correlation matrix must be 1!")
    if (any(c(err.y, err.x)<0)) stop("Fatal error: uncertainties err.y and err.x cannot be negative!")
    
    err.z <- c(err.y, err.x)                      # combined error standard deviation vector
    rzz <- cbind(c(1,rxy), rbind(t(rxy), rxx))    # combine correlations to build correlation matrix of Z
    Szz <- (err.z %*% t(err.z)) * rzz             # covariance matrix of Z
  }
  
  # and unpack Szz into Sxx and Sxy
  Sxy <- as.vector(Szz[,1])                       # covariances of X and y are in the first column of Szz...
  Sxy <- Sxy[-1]                                  # ... with the first element removed (that is err.y^2)
  Sxx <- Szz[-1, -1]                              # covariance matrix of x's is Szz minus the first row and column
  
  
  #############################################################
  # now find lambda as the smallest root of det(mZZ-lambda*Szz)
  
  # First, if there are any constants in X, we need to replace the resulting rows and columns of
  # zeroes in mZZ with 1's, and the corresponding rows and columns of Szz with zeroes.  Otherwise 
  # the determinant in find_root will be badly distorted.
  
  if (sum(diag(mZZ)==0)>0) {
    constants <- which(diag(mZZ)==0)    # these X's are constants
    mZZ <- mZZ[ ,-constants]
    mZZ <- mZZ[-constants, ]
    Szz <- Szz[ ,-constants]
    Szz <- Szz[-constants, ]
  }
  
  
  lambda <- find_root(mZZ, Szz)         # now call find_root to find the smallest root
  
  
  ##################################################
  # now do Fuller's corrections to the crossproducts
  if (uncert=="relative") {                  # Fuller's eq. 2.5.14  
    adj_MXX <- MXX - lambda*(1-alpha/n)*Sxx 
    adj_MXY <- MXY - lambda*(1-alpha/n)*Sxy
  } else {                                   # Fuller's eq. 2.5.13  
    adj_MXX <- MXX - (lambda-alpha/n)*Sxx 
    adj_MXY <- MXY - (lambda-alpha/n)*Sxy
  }
  
  
  #############################
  # now solve the linear system
  inv_adj_MXX <- solve(adj_MXX)                 # invert the crossproduct matrix of the x's
  beta <- inv_adj_MXX %*% adj_MXY               # multiply by the xy crossproduct vector to get the regression slopes
  e    <- y - (X %*% beta)                      # residuals 
  s2e  <- as.vector(sum(e^2) / (n-k))           # residual variance, corrected for degrees of freedom
  beta_cov <- s2e / n * inv_adj_MXX             # covariance matrix of beta: note that we need the factor of 1/n because MXX is a mean square rather than a sum of squares
  se <- sqrt(diag(beta_cov))                    # standard errors of beta coefficients
  
  
  #######################
  # assign variable names
  names(beta) <- xnames            
  names(se) <- xnames
  colnames(beta_cov) <- xnames
  rownames(beta_cov) <- xnames
  
  
  #########################
  # shape and return output
  return(list(
    beta=beta,                  # regression coefficients, corrected for errors in explanatory variables
    se=se,                      # standard errors
    s2e=s2e,                    # residual variance
    beta_cov=beta_cov,          # covariance matrix
    e=e,                        # residuals
    lambda=lambda               # lambda parameter from Fuller's procedure
  ))
  
} 
###########
# end eiv #
###########






###################################################
# ols: performs ordinary least squares regression #
###################################################
ols <- function(
    y,                                  # vector of response variable
    X,                                  # matrix of one or more explanatory variables
    no.intercept = FALSE                # flag for whether regression should be without intercept
)
  
{ 
  
  # receives as input:
  #
  #     y             vector containing the response variable (n rows, 1 column)
  #
  #     X             matrix, data table, or data frame containing one or more explanatory variables (n rows, k columns)
  #                       Column names, if available, are retained for output; otherwise names "x1", "x2", etc. are assigned.
  #
  #                   y and X must have the same number of rows.  NA values are allowed, but only complete cases will be analyzed.
  #
  #     no.intercept  if TRUE, no intercept is fitted.  Removing the intercept can help in constraining the coefficient of (maximum 1!)
  #                       explanatory variable that is constant or nearly constant.
  
  
  # returns a list as output, with the following elements:
  #
  #     beta          OLS regression coefficients
  #
  #     se            standard errors of regression coefficients
  #
  #     s2e           residual variance of regression 
  #
  #     beta_cov      covariance matrix of regression coefficients
  #
  #     e             vector of residuals
  
  
  
  
  ##############
  # check inputs
  
  y <- as.matrix(y)            # make these matrices in case they were provided as vectors, data frames, etc.
  X <- as.matrix(X)
  
  if (NROW(y)!=NROW(X)) stop("fatal error in ols: X and y have different n")
  
  k <- NCOL(X)                  # number of X variables in our problem
  
  if(is.null(colnames(X))) xnames <- paste0("x",1:k)   # copy names of x variables if they are given; otherwise name them x1, x2, etc.
  else xnames <- colnames(X)
  
  # now cbind X and y to make Z, and retain only complete cases
  Z <- cbind(y, X)
  y <- y[complete.cases(Z),]
  X <- X[complete.cases(Z),]
  Z <- Z[complete.cases(Z),]
  
  n <- sum(complete.cases(Z))
  if (n<(k+1)) stop("Fatal error in ols: too few complete cases -- no degrees of freedom left!")
  
  
  #############################################
  # add dummy variable for intercept if desired
  
  if (no.intercept==FALSE) {
    k <- k+1                                 # the width of our design matrix will now expand by 1
    X <- as.matrix(cbind(X,1))               # add a column of 1's to the design matrix X
    xnames <- c(xnames, "intercept")         # and lastly, add "intercept" to the column names
  }
  
  
  ###############################
  # calculate mean crossproducts
  MXX <- crossprod(X)/n
  MXY <- crossprod(X,y)/n
  
  
  #############################
  # now solve the linear system
  inv_MXX <- solve(MXX)                         # invert the covariance matrix of the x's
  beta <- inv_MXX %*% MXY                       # multiply by the xy covariance vector to get the regression slopes
  e    <- y - (X %*% beta)                      # residuals
  s2e  <- as.vector(sum(e^2)) / (n-k)           # residual mean square, corrected for degrees of freedom
  beta_cov <- s2e / n * inv_MXX                 # covariance matrix of beta: note that we need the factor of 1/n because MXX is a mean square rather than a sum of squares
  se <- sqrt(diag(beta_cov))                    # standard errors of beta coefficients
  
  
  
  #######################
  # assign variable names
  names(beta) <- xnames            
  names(se) <- xnames
  colnames(beta_cov) <- xnames
  rownames(beta_cov) <- xnames
  
  
  #########################
  # shape and return output
  return(list(
    beta=beta,                  # regression coefficients, corrected for errors in explanatory variables
    se=se,                      # standard errors
    s2e=s2e,                    # residual variance
    beta_cov=beta_cov,          # covariance matrix
    e=e                         # residuals
  ))
  
} 
###########
# end ols #
###########











##########################################################################################
# eemma: performs ensemble end-member mixing analysis, with propagation of uncertainties #
##########################################################################################
eemma <- function(
    y,                                  # vector containing tracer time series in mixture
    X,                                  # vector or matrix with tracer time series in one or more end-members
    complete = TRUE,                    # flag for whether set of end-members is complete (TRUE), or whether there are unmeasured end-members (FALSE)
    memory = FALSE,                     # flag for whether eemma should quantify the fraction of the mixture in one time step that is carried over from the previous time step
    no.intercept = complete,            # should the intercept be omitted?  TRUE by default (and unnecessary) if complete==TRUE, and FALSE otherwise
    method = c("OLS", "EIV"),           # should solution be by ordinary least squares (OLS) or the error-in-variables method (EIV)?
    err.all = 1,                        # measurement error in X and y, expressed as standard deviation
    err.y = err.all,                    # measurement error in y, expressed as standard deviation
    err.x = rep(err.all, NCOL(X)),      # vector of measurement errors in x's, expressed as standard deviation
    sercorr.err.y = 0,                  # serial correlation in *errors* in y (note: NOT serial correlation in y values themselves)
    rxx = 0,                            # correlation matrix of measurement errors in x's (correlations, not covariances!)
    rxy = 0,                            # correlation between measurement errors of errors in y and each of the x's (correlations, not covariances!)
    uncert = c("relative", "absolute"), # are uncertainties absolute, or relative?
    alpha = 1                           # parameter with value between 1 and 2, inclusive
)

{ 
  
  # receives as input:
  #
  #     y             vector containing sequence of tracer measurements in mixture (n rows, 1 column)
  #
  #     X             matrix, data table, or data frame containing sequence of tracer measurements in one or more end-members (n rows, k columns)
  #                       Column names, if available, are retained for output; otherwise names "em1", "em2", etc. are assigned.
  #
  #                   y and X must have the same number of rows.  NA values are allowed, but only complete cases will be analyzed.
  #
  #     complete      flag for whether the set of end-members is complete (TRUE), in which case their mixing fractions should sum to 1, or 
  #                       whether there are unmeasured end-members (FALSE), in which case the mixing fractions are not constrained to sum to 1,
  #                       and the sum of the unmeasured end-members is estimated from the residuals
  #
  #     memory        flag for whether eemma should estimate a first-order memory coefficient that quantifies how much of the mixture in one time step
  #                       is carried over from the previous time step.  FALSE by default.
  #
  #     no.intercept  if TRUE, no intercept is fitted.  Removing the intercept can help in constraining the coefficients of explanatory
  #                       variables that are constant or nearly constant, if their means are substantially different.  When complete==TRUE,
  #                       intercept is automatically omitted and no.intercept has no effect.
  # 
  #     method        should mixing fractions be estimated using ordinary least squares (OLS, the default), or by the methods of moments approach to 
  #                       errors in variables (EIV)?  OLS estimates of mixing fractions will be better constrained and yield smaller 
  #                       residuals, but are potentially vulnerable to biases from measurement and sampling errors in the end-member tracer
  #                       concentrations (the "errors in variables" problem).  Where measurement and sampling errors are a substantial
  #                       fraction of the total variance of the explanatory variables (the end-member tracer data), the EIV method will 
  #                       typically yield mixing fractions that have much smaller biases, but with somewhat greater uncertainty.  Can be abbreviated.
  #
  #                       The EIV method requires estimates of the measurement/sampling errors in the end-members and mixture.  
  #                       If no estimates are provided, all measurement/sampling uncertainties are assumed to be equal and uncorrelated, 
  #                       with the magnitude of these uncertainties estimated within the EIV procedure itself.
  #
  #     None of the following parameters need to be set for OLS regression.  They are only used for taking account of 
  #     measurement/sampling errors in the end-members.
  #
  #     err.all       standard deviation of measurement/sampling errors in x's and y.  If only err.all is provided, it is assumed to 
  #                       apply to both y and all the x's.  Default value is 1, with (by default) implicit rescaling as invoked by 
  #                       the default uncert=="relative" flag, below.
  #
  #     err.y         standard deviation of measurement/sampling error in y.  If provided, supersedes err.all.
  #
  #     err.x         standard deviation(s) of measurement/sampling errors in x's.  If provided, supersedes err.all.
  #                       Must be a vector of length NCOL(X), or a single value.
  #                       If a single value rather than a vector is supplied, it will be copied for all of the x's.
  #
  #     sercorr.err.y serial correlation in errors in y.  Must be a single value between -1 and 1.  
  #                       Note that this is not the serial correlation in the values of y themselves.
  #
  #     rxx           correlation matrix of measurement/sampling errors in X (k by k matrix).  
  #                       If a single value is provided, a correlation matrix is constructed using that value for all cross-correlations.
  #                       If a matrix is provided, it must be positive definite with values between -1 and 1.
  #                       Any values !=1 on the diagonal will be replaced with 1's.
  #                       If rxx is not specified, all cross-correlations are assumed to be zero.
  #
  #     rxy           vector of correlations between measurement/sampling errors between y and each of the x's.
  #                       If a single value is provided, a vector is constructed using that value for all cross-correlations.
  #                       If rxy is not specified, these cross-correlations are assumed to be zero.
  #
  #     uncert        flag for whether the relative, or absolute, values of the uncertainties have been specified.  
  #                       If uncert=="relative", uncertainties are expressed as relative values (for example, if uncertainty 
  #                       in the mixture measurements is twice the uncertainty in the end-member measurements, err.y and err.x 
  #                       can be any numbers in the ratio of 2:1) and these are then implicitly rescaled in the fitting procedure.
  #                       Can be abbreviated.  Default is "relative".  
  #
  #     alpha         alpha parameter for Fuller's EIV method.  Must be positive. 
  #                       Lower values give smaller bias but larger parameter variance; higher values give more bias but smaller parameter variance.
  #                       Values between 1 and 4 have performed well in benchmark tests.
  
  
  # returns a list as output, with the following elements:
  #
  #     coeffs        a t-table for the mixing fractions and their standard errors
  #
  #     f             estimated mixing fractions (also for the unknown source "unk", if no.intercept==FALSE)
  #
  #     se            standard errors of mixing fractions (also for the unknown source "unk", if no.intercept==FALSE)
  #
  #     f_cov         covariance matrix of mixing fractions (including the unknown source "unk", if no.intercept==FALSE)
  #
  #     RMSE          root-mean-square error (square root of residual variance)
  #
  #     e             residuals
  #
  #     timeseries    matrix of time series of the mixture and measured end-members (and intercept and unknown source "unk", if no.intercept==FALSE)
  #
  #     lambda        if complete == FALSE, the following are also returned
  
  
  
  
  ############################################
  uncert <- match.arg(uncert)  # uncertainties relative or absolute?
  method <- match.arg(method)  # method OLS or EIV?
  
  
  # correlation matrix of measurement errors in x's (correlations, not covariances!)
  if (length(rxx)==1) {                              # if rxx is supplied as a single value
    rxx <- matrix(rxx, nrow=NCOL(X), ncol=NCOL(X))   # copy this for all the cross-correlations in the rxx matrix
    diag(rxx) <- 1                                   # and set the diagonal to 1
  }

  if (NCOL(rxx)!=NCOL(X)) stop("Fatal error in eemma: rxx has wrong size")
  if (NROW(rxx)!=NCOL(X)) stop("Fatal error in eemma: rxx has wrong size")
  if (any(abs(rxx)>1)) stop("Fatal error in eemma: all values of rxx must be between -1 and 1, inclusive")
  if (!isSymmetric(rxx)) stop("Fatal error in eemma: rxx must be symmetric")
  if (min(eigen(rxx, only.values=TRUE)$values)<=0) stop("Fatal error in eemma: rxx must be positive definite")
  
  
  # correlation between measurement errors of errors in y and each of the x's (correlations, not covariances!)
  if (length(rxy)==1) rxy <- rep(rxy, NCOL(X))
  if (length(rxy)!=NCOL(X)) stop("Fatal error in eemma: rxy has wrong size")
  
  if (any(abs(rxy)>1)) stop("Fatal error in eemma: all values of rxy must be between -1 and 1, inclusive")
  
  if (abs(sercorr.err.y)>=1) stop("Fatal error in eemma: sercorr.err.y must be between -1 and 1")
  
  if (alpha<=0) stop("Fatal error in eemma: alpha must be positive")

  y <- as.matrix(y)            # make these matrices in case they were provided as vectors, data frames, etc.
  X <- as.matrix(X)
  
  if (NROW(y)!=NROW(X)) stop("fatal error in eemm: data series have different n")
  
  k <- NCOL(X)                  # number of end-members in our problem
  
  if(is.null(colnames(X))) xnames <- paste0("EM", 1:k)   # copy names of x's if they are given; otherwise name them em1, em2, etc. (as abbreviation for "end-member")
  else xnames <- colnames(X)
  
  if (length(err.x)==1) err.x <- rep(err.x, NCOL(X))    # if err.x was specified as a single value, copy it for all the X variables
  else if (length(err.x)!=NCOL(X)) stop("Fatal error in eemma: err.x must be a single value, or a vector of the same length as NCOL(X)")
  
  
  if (memory==TRUE) {                  # if we are including a first-order memory term
    X <- cbind(X, dplyr::lag(y,1))     # add the lag of y as the right-most column of X  
    xnames <- c(xnames, "lagmix")      # add an appropriate name
    colnames(X) <- xnames
    err.x <- c(err.x, err.y)           # mirror the uncertainty in y here
    k <- k+1                           # increment k
    rxx <- rbind( cbind(rxx, rxy*sercorr.err.y), c(rxy*sercorr.err.y,1) )     # adjust the rxx matrix accordingly as well
    rxy <- c(rxy, sercorr.err.y)       # add serial correlations in *errors* (note: NOT serial correlation in values of y itself) to the rxy vector
  }

  
  
  #############################################################
  # now cbind y and X to make Z, and retain only complete cases
  Z <- cbind(y, X)
  y <- y[complete.cases(Z),]
  X <- X[complete.cases(Z),]
  Z <- Z[complete.cases(Z),]
  
  n <- sum(complete.cases(Z))
  
  
  ####################################################
  # now compile error variance and covariance matrices
  
  if (det(rxx)<=0) stop("Fatal error: error correlation matrix must be positive definite!")
  if (any(diag(rxx)!=1)) stop("Diagonal elements of correlation matrix must be 1!")
  if (any(c(err.x, err.y)<0)) stop("Fatal error: uncertainties err.y and err.x cannot be negative!")
  
  err.z <- c(err.y, err.x)                     # combined error standard deviation vector
  
  rzz <- cbind(c(1,rxy), rbind(t(rxy), rxx))   # combine correlations to build correlation matrix of Z
  if (min(eigen(rzz, only.values=TRUE)$values)<=0) stop("Fatal error in eemma: rzz must be positive definite")
  Szz <- (err.z %*% t(err.z)) * rzz            # covariance matrix of Z
  
  if (min(eigen(Szz, only.values=TRUE)$values)<=0) stop("Fatal error in eemma: Szz must be positive definite")
  
  
  ########################################################################################################################
  # if complete==TRUE, subtract one column from X and y, transform Z and Szz, run analysis,then back-transform the results
  if (complete==TRUE) {
    
    if (no.intercept==FALSE) warning("Warning: eemma option no.intercept==FALSE is being ignored because it is inconsistent with complete==TRUE")
    
    if (NCOL(X)==1) stop("Fatal error: complete==TRUE requires two or more end-members, or one end-member if memory==TRUE!")
    
    sds <- apply(X, 2, sd)                     # calculate std deviation of each column of X
    if (sum(sds==0)>2) stop("Fatal error: if complete==TRUE, there cannot be more than two constant end-members")
    
    
    ###########################################################
    # now we're going to transform y and X all at the same time.  We subtract the last column of X from everything else (mathematically it makes no difference which column we pick!)
    tr <- diag(NCOL(Z))              # build transformation matrix starting with identity matrix
    tr[(k+1), ] <- -1                # make row k of transformation matrix -1 (note this is row k+1 of Z matrix)
    tr <- tr[ ,-(k+1)]               # and remove corresponding column of transformation matrix
    
    trz <- Z %*% tr                  # transform the Z matrix
    trSzz <- t(tr) %*% (Szz %*% tr)  # and transform the covariance matrix
    
    
    ###########################################################################
    # now run the regression with the *transformed* y, X, and covariance matrix
    if (method=="EIV") q <- eiv(y=trz[,1], X=trz[,-1], no.intercept=TRUE, Szz=trSzz, uncert=uncert, alpha=alpha)
    else q <- ols(y=trz[,1], X=trz[,-1], no.intercept=TRUE)
    
    
    ###############################################
    # then we back-transform the betas and beta_cov  
    
    if (memory==TRUE) {                         # if our last column was the lag of y (to account for memory effects), we need to transform all the beta's
      
      sumbeta <- sum(q$beta)                    # sum of beta's
      f <- q$beta/sumbeta                       # mixing fractions
      f <- c(f, 1-sumbeta)                      # f[k+1] is not technically a mixing fraction, but instead the memory fraction

      J <- matrix(-q$beta/sumbeta^2, nrow=k-1, ncol=k-1)   # Jacobian matrix for this problem (step 1)
      J <- J + diag(1/sumbeta, nrow=k-1)                   # Jacobian matrix for this problem (step 2)
      J <- rbind(J, rep(-1, n=k-1))                        # Jacobian matrix for this problem (step 3)
      
      f_cov <- J %*% q$beta_cov %*% t(J)        # propagate covariance matrix
      se <- sqrt(diag(f_cov))                   # standard errors of f's
      
    } else {                                    # if we have a conventional "complete" case, without memory
      
      btr <- tr[-1,-1]                          # back-transformation for coefficients omits first row and column (which corresponded to the response variable)
      f <- btr %*% q$beta                       # back-transform beta to get mixing fractions f
      f[k] <- 1+f[k]                            # need to add 1 to the (negative) sum of f's for the end-member that was subtracted
      
      f_cov <- btr %*% q$beta_cov %*% t(btr)    # propagate uncertainties in beta to uncertainties in f's
      se <- sqrt(diag(f_cov))                   # standard errors of f's
      
    } # end if memory==TRUE
    
    
    
    timeseries <- as.matrix(cbind(y, X))
    colnames(timeseries) <- c("mixture", xnames)
    
    
    #########################
    # shape and return output
    names(se) <- xnames
    names(f) <- xnames
    colnames(f_cov) <- xnames
    rownames(f_cov) <- xnames
    
    coeffs = data.frame(f=f, se=se, t_statistic=abs(f/se), p_value=pt(abs(f/se), df=n-k, lower.tail=FALSE))
    
    results <- list(
      coeffs=coeffs,
      f_cov = f_cov,
      RMSE = sqrt(q$s2e),
      f = f,
      se = se,
      e = q$e,
      timeseries=timeseries,
      lambda=NA
    )
    
    if (method=="EIV") results$lambda=q$lambda
    
    return(results)
    
    
  } #end if complete==TRUE  (we do not need an "else" here, because any complete==TRUE cases have already returned and exited)
  
  
  
  
  ###########################################################################################################################
  # if complete==FALSE, we do not not transform the input data, except for centering around the mixture mean ybar (if no.intercept==TRUE)
  
  sds <- apply(as.matrix(X), 2, sd)  # calculate std deviation of each column of X
  if (no.intercept) {
    if (sum(sds==0)>1) stop("Fatal error: if complete==FALSE, and no.intercept==TRUE, there can be a maximum of one constant end-member")
  } else {
    if (sum(sds==0)>0) stop("Fatal error: if complete==FALSE and no.intercept==FALSE, there can be no constant end-members")
  }
  
  
  # If we have no intercept, subtract the average mixture concentration from y and all of the x's.  This "centers" the problem around the average mixture, avoiding
  # artifacts that can arise when the measurement scale puts all of the values far from zero, and an intercept is not used.
  if (no.intercept==TRUE) ybar <- mean(y) else ybar <- 0
  y_ybar <- y-ybar
  X_ybar <- X-ybar
  
  
  #######################################################
  # now run the regression with the *transformed* y and X (centering is just a translation, without rescaling, so Szz is not affected)
  if (method=="EIV") q <- eiv(y=y_ybar, X=X_ybar, no.intercept=no.intercept, Szz=Szz, uncert=uncert, alpha=alpha)
  else q <- ols(y=y_ybar, X=X_ybar, no.intercept=no.intercept)
  
  
  #############################
  # now transform outputs
  
  if (memory==TRUE) {  # if we have included the lag of y to estimate the first-order memory effect
    
    
    if (no.intercept==TRUE) {  #if we have no intercept, then beta[k] is the "memory" term
      
      f <- q$beta/(1-q$beta[k])      # mixing fractions
      f[k] <- q$beta[k]              # not really a mixing fraction; this is the "memory" coefficient
      
      J <- diag(1/(1-q$beta[k]), nrow=k-1)       # Jacobian matrix (step 1)
      J <- cbind(J, q$beta[-k]/(1-q$beta[k])^2)  # Jacobian matrix (step 2)
      J <- rbind(J, c(rep(0, times=k-1),1))      # Jacobian matrix (step 3)
      
      f_cov <- J %*% q$beta_cov %*% t(J)            # propagation of covariance matrix
      se <- sqrt(diag(f_cov))        # standard errors
      e <- q$e                       # residuals
      
      timeseries <- as.matrix(cbind(y, X))
      colnames(timeseries) <- c("mixture", xnames)
      
    } else {          # if no.intercept==FALSE, calculate possible unknown source and its associated standard errors
      
      # here, q$beta[k] is the memory term and q$beta[k+1] is the intercept
      
      f <- rep(0, times=k+2)
      f[1:(k-1)] <- q$beta[1:(k-1)]/(1-q$beta[k])      # mixing fractions
      f[k] <- q$beta[k]                                # not really a mixing fraction; this is the "memory" coefficient
      f[k+1] <- q$beta[k+1]                            # also not really a mixing fraction; this is the intercept
      f[k+2] <- 1-sum(q$beta[1:(k-1)])/(1-q$beta[k])   # estimate of mixing fraction of unknown source
      
      J <- diag(1/(1-q$beta[k]), nrow=k-1)       # Jacobian matrix (step 1)
      J <- cbind(J, q$beta[1:(k-1)]/(1-q$beta[k])^2)  # Jacobian matrix (step 2) column for effects of memory coefficient
      J <- rbind(J, c(rep(0, times=k-1),1))      # Jacobian matrix (step 3) row for propagation of uncertainty in memory coefficient (nothing happens here!)
      J <- cbind(J, rep(0, times=k))             # Jacobian matrix (step 4) column for effects of intercept (there aren't any)
      J <- rbind(J, c(rep(0, times=k),1))        # Jacobian matrix (step 5) row for propagation of uncertainty in intercept (nothing happens here!)
      J <- rbind(J, c(rep(-1/(1-q$beta[k]), times=k-1), sum(q$beta[1:(k-1)])/(1-q$beta[k])^2, 0))    # Jacobian matrix (step 6) new row for propagation of uncertainty in f_unk
      
      f_cov <- J %*% q$beta_cov %*% t(J)                    # propagate covariance matrix
      
      se <- sqrt(diag(f_cov))                               # standard error including unknown source (also accounting for parameter covariances)

      xnames <- c(xnames, "intercept", "unk")               # add labels for an additional parameter (intercept) and an estimate of the mixing fraction of an unknown end-member

      e <- q$e                                              # residuals
      
      unk <- ybar + (f[k+1] + e)/(f[k+2]*(1-f[k]))          # time series of unknown source
      
      du.dint <- 1/(f[k+2]*(1-f[k]))                        # derivative of unk w.r.t. intercept (f[k+1])  - for standard error calculation
      du.dfu <- -(f[k+1] + e)/(f[k+2]*(1-f[k]))^2*(1-f[k])  # derivative of unk w.r.t. f_u (f[k+2])
      du.dfm <- (f[k+1] + e)/(f[k+2]*(1-f[k]))^2*f[k+2]     # derivative of unk w.r.t. f_memory (f[k])
            
      # standard error of unknown source time series (1st-order 2nd-moment error propagation)
      unk_se <- sqrt( (se[k]*du.dfm)^2 + (se[k+1]*du.dint)^2 + (se[k+2]*du.dfu)^2 + 2*du.dfm*du.dint*f_cov[k, k+1] + 2*du.dfm*du.dfu*f_cov[k, k+2] + 2*du.dint*du.dfu*f_cov[k+1, k+2])
      
      timeseries <- as.matrix(cbind(y, X, f[k+1]+ybar, unk, unk_se))
      colnames(timeseries) <- c("mixture", xnames, "unk_std_err")
      
    } # end if no.intercept == FALSE

        
  } else {  # if we aren't estimating memory effects
    
    f <- q$beta                    # mixing fractions
    f_cov <- q$beta_cov            # covariance matrix
    se <- sqrt(diag(f_cov))        # standard errors
    e <- q$e                       # residuals
    
    
    ##############################################################################################
    # if no.intercept==FALSE, calculate possible unknown source and its associated standard errors
    if (no.intercept==FALSE) {
      xnames <- c(xnames, "intercept", "unk")               # if we estimated an intercept, we have an additional parameter (intercept) and an estimate of the mixing fraction of an unknown end-member
      
      f <- c(f, 1-sum(f[1:k]))                              # the unknown source has a mixing fraction of 1-sum(f's) (not including the intercept)
      
      f_tr <- rbind(diag(k+1), c(rep(-1,k), 0))             # transformation matrix from f's without unknown source, to f's with unknown source
      
      f_cov <- f_tr %*% (f_cov %*% t(f_tr))                 # updating f_cov to include f of unknown source
      
      se <- sqrt(diag(f_cov))                               # standard error including unknown source (also accounting for parameter covariances)
      
      unk <- ybar + (f[k+1] + e)/f[k+2]                     # time series of unknown source
      
      # standard error of unknown source time series (1st-order 2nd-moment error propagation)
      unk_se <- sqrt((se[k+1]/f[k+2])^2 + (se[k+2]*(f[k+1] + e)/f[k+2]^2)^2 - 2*f_cov[k+1,k+2]*(f[k+1] + e)/f[k+2]^3)

      timeseries <- as.data.frame(cbind(y, X, f[k+1]+ybar, unk, unk_se))
      colnames(timeseries) <- c("mixture", xnames, "unk_std_err")
    } else {
      timeseries <- as.data.frame(cbind(y, X))
      colnames(timeseries) <- c("mixture", xnames)
    }
    
  }
  
  
  #########################
  # shape and return output
  names(se) <- xnames
  names(f) <- xnames
  colnames(f_cov) <- xnames
  rownames(f_cov) <- xnames
  
  # compute and compile t-table.  Note that p-value is *not* multiplied by 2 because these are one-tailed p's
  # for tests with null hypothesis of f<=0 (not f=0, which would be a conventional two-tailed test).
  coeffs = data.frame(f=f, se=se, t_statistic=abs(f/se), p_value=pt(abs(f/se), df=n-k, lower.tail=FALSE))
  
  results <- list(
    coeffs=coeffs,
    f_cov = f_cov,
    RMSE = sqrt(q$s2e),
    f = f,
    se = se,
    e = q$e,
    timeseries=timeseries,
    lambda=NA
  )
  
  
  if (method=="EIV") results$lambda=q$lambda
  
  return(results)
  
} 
#############
# end eemma #
#############


