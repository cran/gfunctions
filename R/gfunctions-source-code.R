###############################################################################
##
## Contents:
##
## 1 gsummary()
##   gsummary.data.frame()
##   gsummary.lm()
##   gsummary.glm()
## 2 glag()
## 3 ghist(): for the future
##
###############################################################################


###############################################################################
## 1 gsummary()
###############################################################################

##==================================================================
## create generic:
gsummary <- function(object, ...){ UseMethod("gsummary") }

##==================================================================
## gsummary.default():
gsummary.default <- function(object, ...)
{
  summary(object, ...)
}

##==================================================================
## gsummary.data.frame():
##
## The gsummary.data.frame() function summarise the main properties
## of a data frame in a way similar to that of STATA
##  
## Arguments:
##
##    object      An object of class 'data.frame'
##
## Example:
##
##    set.seed(123)
##    x <- rnorm(20); y <- rnorm(20); z <- rnorm(20)
##    mydataframe <- as.data.frame(cbind(x,y,z))
##    gsummary(mydataframe)
##
##==================================================================
##create the 'gsummary.data.frame' S3 method:
gsummary.data.frame <- function(object, ...)
{
  varnames <- colnames(object)
  ncols <- length(varnames)
  result <- matrix(NA, ncols, 5)
  colnames(result) <- c("Obs", "Mean", "Std. Dev.", "Min", "Max") 
  rownames(result) <- varnames  
  
  ##fill:
  for(i in 1:ncols){
    result[i,"Obs"] <- length( which( !is.na(object[,i]) ) )
  }
  result[,"Mean"] <- colMeans(object, na.rm=TRUE)    
  result[,"Std. Dev."] <- sqrt(diag(var(object, na.rm=TRUE)))
  for(i in 1:ncols){ result[i,"Min"] <- min(object[,i], na.rm=TRUE) }
  for(i in 1:ncols){ result[i,"Max"] <- max(object[,i], na.rm=TRUE) }

  ##print result:
  print(result)
}

##==================================================================
## gsummary.lm():
##
## The gsummary.lm() function provides a print of the estimation
## results similar to that of STATA for an object of class 'lm'.
##  
## Arguments:
##
##    object        An object of class 'lm'
##
##    vcov.type     A character, either "ordinary", "robust" or "hac".
##                  If "ordinary", the default, then ordinary standard
##                  errors are used. If "robust", then the
##                  heteroscedasticity robust standard errors of White
##                  (1980) are used. If "hac", then the heteroscedasticity
##                  and autocorrelation consistent standard errors of
##                  Newey and West (1987) are used.
##            
##    confint.level A number between 0 and 1, the confidence level for the
##                  confidence intervals of the coefficient estimates. If
##                  NULL, then confidence intervals are not computed.
##
## Example:
##
##    set.seed(123)
##    z <- rnorm(20); u <- rnorm(20); y <- rnorm(20)
##    mymodel <- lm(y ~ z + u)
##    gsummary(mymodel)
##    gsummary(mymodel, vcov.type="robust")
##    gsummary(mymodel, confint.level=0.90)
##
##==================================================================
##create the 'gsummary.lm' S3 method:
gsummary.lm <- function(object, vcov.type = c("ordinary", "robust", "hac"),
  confint.level = 0.95, ...)
{
  ##record names
  ##------------
  x <- object
  xName <- deparse(substitute(x))
  yName <- names(x$model)[1]

  ##check if x is 'lm' object
  ##-------------------------
  if( !is(x,"lm") ){
    warning(paste0("'", xName, "' not of class 'lm'"))
  }
    
  ##vcov.type argument
  ##------------------
  types <- c("ordinary", "robust", "hac")
  whichType <- charmatch(vcov.type[1], types)
  vcov.type <- types[ whichType ]
  #Or: vcov.type <- match.arg(vcov.type)

  ##check 'confint.level' argument
  ##------------------------------
  if( !is.null(confint.level) ){
    if( confint.level <= 0 || confint.level >= 1 ) stop("'confint.level' must be between 0 and 1")
  }
  
  ##record stuff
  ##------------
  coefs <- coef(x)
  vcovmat <- vcov(x)
  resids <- resid(x)
  n <- nobs(x)
  k <- length(coefs)  

  ##x has intercept?
  ##----------------
  if( length(coefs)==0 ){
    xHasIntercept <- FALSE
  }else{
    tmp <- model.matrix(x)
    attr(tmp, "assign") <- NULL #remove 'assign' attribute
    xHasIntercept <- ifelse(colnames(tmp)[1] == "(Intercept)", TRUE, FALSE)
  }
                                
  ##print header
  ##------------
  cat("\n")
  cat("Date:", date(), "\n")
  #cat("Name:", xName, "\n")
  cat("Dependent var.:", yName, "\n")
  cat("Coefficient covariance (vcov):", vcov.type, "\n")
  cat("Number of obs. (nobs):", n, "\n")  

  ##goodness-of-fit
  ##---------------

  ##create:
  GOFresults <- matrix(NA, 3, 8)
  rownames(GOFresults) <- c("Model: ", "Residual: ", "Total: ")
  
  ##fill SS part:
  RSS <- sum(residuals(x)^2)
  y <- as.vector(x$model[,1])
  ymean <- mean(y)
  TSS <- sum( (y-ymean)^2 )
  ESS <- TSS-RSS
  GOFresults[1,1] <- ESS #fill column 1
  GOFresults[2,1] <- RSS
  GOFresults[3,1] <- TSS
  GOFresults[1,2] <- k - as.numeric(xHasIntercept) #fill column 2
  GOFresults[2,2] <- n-k
  GOFresults[3,2] <- n - as.numeric(xHasIntercept)

  ##fill gof part:
  R2 <- 1 - RSS/TSS
  GOFresults[1,8] <- R2
  GOFresults[2,8] <- 1 - ( (1-R2)*(n-1)/(n-k) )
  GOFresults[3,8] <- sqrt( RSS/(n-k) )

  ##convert to data frame, complete
  GOFresults <- as.data.frame(GOFresults)
  colnames(GOFresults) <- c("SS", "df", "", "", "", "", "", "")
  GOFresults[,7] <- c("R-squared =", "Adj. R-squared =", "Root MSE =")
  GOFresults[1:3,3:6] <- ""
        
  ##print:
  cat("\n")
  cat("Goodness-of-fit:", "\n")
  cat("\n")
  print(GOFresults)
  #printCoefmat(GOFresults, signif.stars=FALSE)

  ##F-statistic
  Fresults <- matrix(NA, 1, 5)
  noOfSlopeCoefs <- k - as.numeric(xHasIntercept)
  if( noOfSlopeCoefs > 0 ){
    Fstatistic <- ((TSS-RSS)/noOfSlopeCoefs)/(RSS/(n-k))
    Fresults[1,1] <- Fstatistic
    Fresults[1,5] <- pf(Fstatistic, noOfSlopeCoefs, n-k, lower.tail=FALSE)
  }
  Fresults <- as.data.frame(Fresults)
  Fresults[1,4] <- c("Prob > F =")
  Fresults[1,2:3] <- ""
  colnames(Fresults) <- c("", "", "", "", "")
  Fstatname <- paste0("F(", noOfSlopeCoefs, ", ", n-k, ") =")
  rownames(Fresults) <- Fstatname
  #rownames(Fresults) <- c("F-statistic =")

  ##print:
  print(Fresults)
    
  ##main table:
  ##-----------
  
  ##create:
  xnames <- c("coefs", "std.error", "t-value", "p-value")
  xresults <- matrix(NA, k, length(xnames))
  colnames(xresults) <- xnames
  rownames(xresults) <- names(coefs)
  
  ##vcov matrix:
  if( vcov.type=="ordinary" ){ vcovMat <- vcov(x) }
  if( vcov.type=="robust" ){ vcovMat <- vcovHC(x, type="HC") }
  if( vcov.type=="hac" ){ vcovMat <- NeweyWest(x) }  
  
  ##fill:
  xresults[,"coefs"] <- coef(x)
  stderrors <- sqrt(diag(vcovMat))
  xresults[,"std.error"] <- stderrors
  tstats <- xresults[,"coefs"]/xresults[,"std.error"]
  xresults[,"t-value"] <- tstats
  xresults[,"p-value"] <- pt(abs(tstats), n-k, lower.tail=FALSE)*2
    
  ##confidence intervals:
  if( !is.null(confint.level) ){
    confintMat <- matrix(NA, length(coefs), 2)
    ci <- format(100*confint.level, nsmall=0)
    lowerName <- paste0("lower ", ci, "%")
    upperName <- paste0("upper ", ci, "%")
    colnames(confintMat) <- c(lowerName, upperName)
    halfalpha <- (1 - confint.level)/2
    tcritval <- qt(1-halfalpha, n-k)
    uppervals <- as.numeric(coefs + tcritval*stderrors)
    lowervals <- as.numeric(coefs - tcritval*stderrors)
    confintMat[,1] <- lowervals
    confintMat[,2] <- uppervals
    rownames(confintMat) <- rownames(xresults)
    confintMat <- as.data.frame(confintMat)
    xresults <- cbind(xresults, confintMat)
  }
  
  ##print main table:
  cat("\n")
  cat("Estimation results:", "\n")
  cat("\n")
  print(xresults)
  #printCoefmat(xresults, signif.stars=signif.stars)

} #close gsummary.lm()


##==================================================================
## gsummary.glm():
##
## The gsummary.glm() function provides a print of the estimation
## results similar to that of STATA for an object of class 'glm'.
##  
## Arguments:
##
##    object        An object of class 'lm'
##
##    vcov.type     A character, either "ordinary", "robust" or "hac".
##                  If "ordinary", the default, then ordinary standard
##                  errors are used. If "robust", then the
##                  heteroscedasticity robust standard errors of White
##                  (1980) are used. If "hac", then the heteroscedasticity
##                  and autocorrelation consistent standard errors of
##                  Newey and West (1987) are used.
##            
##    confint.level A number between 0 and 1, the confidence level for the
##                  confidence intervals of the coefficient estimates. If
##                  NULL, then confidence intervals are not computed.
##
## Example:
##
##    set.seed(123)
##    y <- as.numeric( c(runif(20)-0.5) > 0 ); x <- rnorm(20); z <- rnorm(20)
##    mymodel <- glm(y ~ x + z, family=binomial)
##    gsummary(mymodel)
##    gsummary(mymodel, vcov.type="robust")
##    gsummary(mymodel, confint.level=0.90)
##
##==================================================================
##create the 'gsummary.glm' S3 method:
gsummary.glm <- function(object, confint.level = 0.95, ...)
{
  ##record names
  ##------------
  x <- object
  xName <- deparse(substitute(x))
  yName <- names(x$model)[1]

  ##check if x is 'glm' object
  ##-------------------------
  if( !is(x,"glm") ){
    warning(paste0("'", xName, "' not of class 'glm'"))
  }
    
#  ##vcov.type argument
#  ##------------------
#  types <- c("ordinary", "robust", "hac")
#  whichType <- charmatch(vcov.type[1], types)
#  vcov.type <- types[ whichType ]
#  #Or: vcov.type <- match.arg(vcov.type)

  ##check 'confint.level' argument
  ##------------------------------
  if( !is.null(confint.level) ){
    if( confint.level <= 0 || confint.level >= 1 ) stop("'confint.level' must be between 0 and 1")
  }
  
  ##record stuff
  ##------------
  coefs <- coef(x)
  vcovmat <- vcov(x)
  resids <- resid(x)
  n <- nobs(x)
  k <- length(coefs)  

  ##x has intercept?
  ##----------------
  if( length(coefs)==0 ){
    xHasIntercept <- FALSE
  }else{
    tmp <- model.matrix(x)
    attr(tmp, "assign") <- NULL #remove 'assign' attribute
    xHasIntercept <- ifelse(colnames(tmp)[1] == "(Intercept)", TRUE, FALSE)
  }
                                
  ##print header
  ##------------
  cat("\n")
  cat("Date:", date(), "\n")
  #cat("Name:", xName, "\n")
  cat("Dependent var.:", yName, "\n")
#  cat("vcov:", vcov.type, "\n")
  cat("Number of obs.:", n, "\n")  

#  ##goodness-of-fit
#  ##---------------
#
#  ##create:
#  GOFresults <- matrix(NA, 3, 8)
#  rownames(GOFresults) <- c("Model: ", "Residual: ", "Total: ")
#  
#  ##fill SS part:
#  RSS <- sum(residuals(x)^2)
#  y <- as.vector(x$model[,1])
#  ymean <- mean(y)
#  TSS <- sum( (y-ymean)^2 )
#  ESS <- TSS-RSS
#  GOFresults[1,1] <- ESS #fill column 1
#  GOFresults[2,1] <- RSS
#  GOFresults[3,1] <- TSS
#  GOFresults[1,2] <- k - as.numeric(xHasIntercept) #fill column 2
#  GOFresults[2,2] <- n-k
#  GOFresults[3,2] <- n - as.numeric(xHasIntercept)
#
#  ##fill gof part:
#  R2 <- 1 - RSS/TSS
#  GOFresults[1,8] <- R2
#  GOFresults[2,8] <- 1 - ( (1-R2)*(n-1)/(n-k) )
#  GOFresults[3,8] <- sqrt( RSS/(n-k) )
#
#  ##convert to data frame, complete
#  GOFresults <- as.data.frame(GOFresults)
#  colnames(GOFresults) <- c("SS", "df", "", "", "", "", "", "")
#  GOFresults[,7] <- c("R-squared =", "Adj.R-squared =", "Root MSE =")
#  GOFresults[1:3,3:6] <- ""
#        
#  ##print:
#  cat("\n")
#  cat("Goodness-of-fit:", "\n")
#  cat("\n")
#  print(GOFresults)
#  #printCoefmat(GOFresults, signif.stars=FALSE)
#
#  ##F-statistic
#  Fresults <- matrix(NA, 1, 5)
#  noOfSlopeCoefs <- k - as.numeric(xHasIntercept)
#  if( noOfSlopeCoefs > 0 ){
#    Fstatistic <- ((TSS-RSS)/noOfSlopeCoefs)/(RSS/(n-k))
#    Fresults[1,1] <- Fstatistic
#    Fresults[1,5] <- pf(Fstatistic, noOfSlopeCoefs, n-k, lower.tail=FALSE)
#  }
#  Fresults <- as.data.frame(Fresults)
#  Fresults[1,4] <- c("Prob > F =")
#  Fresults[1,2:3] <- ""
#  colnames(Fresults) <- c("", "", "", "", "")
#  Fstatname <- paste0("F(", noOfSlopeCoefs, ", ", n-k, ") =")
#  rownames(Fresults) <- Fstatname
#  #rownames(Fresults) <- c("F-statistic =")
#
#  ##print:
#  print(Fresults)
    
  ##main table:
  ##-----------
  
  ##create:
  xnames <- c("coefs", "std.error", "t-value", "p-value")
  xresults <- matrix(NA, k, length(xnames))
  colnames(xresults) <- xnames
  rownames(xresults) <- names(coefs)
  
  ##vcov matrix:
  vcovMat <- vcov(x)
  
  ##fill:
  xresults[,"coefs"] <- coef(x)
  stderrors <- sqrt(diag(vcovMat))
  xresults[,"std.error"] <- stderrors
  tstats <- xresults[,"coefs"]/xresults[,"std.error"]
  xresults[,"t-value"] <- tstats
  xresults[,"p-value"] <- pt(abs(tstats), n-k, lower.tail=FALSE)*2
    
  ##confidence intervals:
  if( !is.null(confint.level) ){
    confintMat <- matrix(NA, length(coefs), 2)
    ci <- format(100*confint.level, nsmall=0)
    lowerName <- paste0("lower ", ci, "%")
    upperName <- paste0("upper ", ci, "%")
    colnames(confintMat) <- c(lowerName, upperName)
    halfalpha <- (1 - confint.level)/2
    tcritval <- qt(1-halfalpha, n-k)
    uppervals <- as.numeric(coefs + tcritval*stderrors)
    lowervals <- as.numeric(coefs - tcritval*stderrors)
    confintMat[,1] <- lowervals
    confintMat[,2] <- uppervals
    rownames(confintMat) <- rownames(xresults)
    confintMat <- as.data.frame(confintMat)
    xresults <- cbind(xresults, confintMat)
  }
  
  ##print main table:
  cat("\n")
  cat("Estimation results:", "\n")
  cat("\n")
  print(xresults)
  #printCoefmat(xresults, signif.stars=signif.stars)

} #close gsummary.glm()


###############################################################################
## 2 glag()
###############################################################################

##==================================================
##create 'glag' generic:
glag <- function(x, ...){ UseMethod("glag") }

##==================================================
##create default S3 method:
glag.default <- function(x, k = 1, pad = TRUE, pad.value = NA, ...)
{
  #check arguments:
  if(k < 1) stop("Lag order k cannot be less than 1")

  #ts and zoo related:
  ists.chk <- is.ts(x)
  iszoo.chk <- is.zoo(x)
  x <- as.zoo(x)
  x.index <- index(x)
  x <- coredata(x)
  isvector <- is.vector(x)
  x <- cbind(x)
  x.n <- NROW(x)
  x.ncol <- NCOL(x)

  #do the lagging:
  x.nmink <- x.n - k
  xlagged <- matrix(x[1:x.nmink,], x.nmink, x.ncol)
  if(pad){
    xlagged <- rbind( matrix(pad.value,k,x.ncol) , xlagged)
  }

  #transform to vector?:
  if(x.ncol==1 && isvector==TRUE){
    xlagged <- as.vector(xlagged)
  }

  #ts or zoo?:
  if(ists.chk || iszoo.chk){
    if(pad){
      xlagged <- zoo(xlagged, order.by=x.index)
    }else{
      xlagged <- zoo(xlagged, order.by=x.index[c(k+1):x.n])
    } #end if(pad)
    if( ists.chk ){ xlagged <- as.ts(xlagged) }
  } #end if(iszoo.chk)

  #out:
  return(xlagged)
} #end glag


###############################################################################
## X ghist()
###############################################################################

###==================================================
###create 'ghist':
#ghist <- function(object, ...){ UseMethod("ghist") }
