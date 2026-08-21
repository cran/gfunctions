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
##    digits      The minimum number of significant digits to be
##                printed (for most numbers)
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
gsummary.data.frame <- function(object, digits=6, ...)
{
  varnames <- colnames(object)
  ncols <- length(varnames)
  
  ##are any of the variables factors?
  isFactor <- rep(NA, ncols)
  for(i in 1:ncols){ isFactor[i] <- is.factor(object[,i]) }

  ##create result object:    
  result <- matrix(NA, ncols, 5)
  colnames(result) <- c("Obs", "Mean", "Std. Dev.", "Min", "Max") 
  rownames(result) <- varnames  
  result <- as.data.frame(result)
  if( any(isFactor) ){
    result <- cbind(result, NA)   
    colnames(result)[6] <- "is.factor()"
    result[,6] <- FALSE
    result[which(isFactor),6] <- TRUE
  }
  
  ##loop over variables:
  for(i in 1:ncols){
    result[i,"Obs"] <- length( which( !is.na(object[,i]) ) )
    if( !isFactor[i] ){
      result[i,"Mean"] <- mean(object[,i], na.rm=TRUE)    
      result[i,"Std. Dev."] <- sqrt(var(object[,i], na.rm=TRUE))
      result[i,"Min"] <- min(object[,i], na.rm=TRUE)
      result[i,"Max"] <- max(object[,i], na.rm=TRUE)
    }
  } #close for(i)

  ##print result:
  print.data.frame(result, digits=digits)
#  format(result, digits=digits, nsmall=digits, scientific=FALSE) 
   
} #close gsummary.data.frame()


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
##    digits        The minimum number of significant digits to be
##                  printed (for most numbers)
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
  confint.level = 0.95, digits=6, ...)
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
  GOFresults <- as.data.frame(GOFresults)
  rownames(GOFresults) <- c("Model: ", "Residual: ", "Total: ")
  colnames(GOFresults) <- c("SS", "df", "", "", "", "", "", "")
  
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
  R2 <- 1 - RSS/TSS #R-squared
  GOFresults[1,8] <- R2
  adjR2 <- 1 - ( (1-R2)*(n-1)/(n-k) ) #adjusted R-squared
  GOFresults[2,8] <- adjR2
  GOFresults[3,8] <- sqrt( RSS/(n-k) )

  ##complete:
  GOFresults[,7] <- c("R-squared =", "Adj. R-squared =", "Root MSE =")
  GOFresults[1:3,3:6] <- ""
        
  ##print:
  cat("\n")
  cat("Goodness-of-fit:", "\n")
  cat("\n")
  print.data.frame(GOFresults, digits=digits)
  #format(GOFresults, digits=digits, nsmall=digits, scientific=FALSE)  

  ##F-statistic
  Fresults <- matrix(NA, 1, 5)
  Fresults <- as.data.frame(Fresults)
  noOfSlopeCoefs <- k - as.numeric(xHasIntercept)
  if( noOfSlopeCoefs > 0 ){
    Fstatistic <- ((TSS-RSS)/noOfSlopeCoefs)/(RSS/(n-k))
    Fresults[1,1] <- format(Fstatistic, digits=digits, scientific=FALSE)
    PvalTxt <- format(pf(Fstatistic, noOfSlopeCoefs, n-k, lower.tail=FALSE),
      digits=digits, scientific=FALSE)
    Fresults[1,5] <- substr(PvalTxt, 1, digits+2)
  }
  Fresults[1,4] <- c("Prob > F =")
  Fresults[1,2:3] <- ""
  colnames(Fresults) <- c("", "", "", "", "")
  Fstatname <- paste0("F(", noOfSlopeCoefs, ", ", n-k, ") =")
  rownames(Fresults) <- Fstatname
  #rownames(Fresults) <- c("F-statistic =")

  ##print:
  print.data.frame(Fresults, digits=digits)
  #format(Fresults, digits=digits, scientific=FALSE)  
  
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
  xresults <- as.data.frame(xresults)
  coefs <- coef(x)
  xresults[,"coefs"] <- format(coefs, digits=digits, scientific=FALSE)
  stderrors <- sqrt(diag(vcovMat))
  xresults[,"std.error"] <- format(stderrors, digits=digits, scientific=FALSE)
  tstats <- coefs/stderrors
  xresults[,"t-value"] <- format(tstats, digits=digits, scientific=FALSE)
  PvalTxt <- format(pt(abs(tstats), n-k, lower.tail=FALSE)*2, digits=digits,
    scientific=FALSE)
  xresults[,"p-value"] <- substr(PvalTxt, 1, digits+2)
#  xresults <- as.data.frame(xresults)
    
  ##confidence intervals:
  if( !is.null(confint.level) ){
    confintMat <- matrix(NA, length(coefs), 2)
    confintMat <- as.data.frame(confintMat)
    ci <- format(100*confint.level, nsmall=0)
    lowerName <- paste0("lower ", ci, "%")
    upperName <- paste0("upper ", ci, "%")
    colnames(confintMat) <- c(lowerName, upperName)
    halfalpha <- (1 - confint.level)/2
    tcritval <- qt(1-halfalpha, n-k)
    uppervals <- as.numeric(coefs + tcritval*stderrors)
    lowervals <- as.numeric(coefs - tcritval*stderrors)
    confintMat[,1] <- format(lowervals, digits=digits, scientific=FALSE)
    confintMat[,2] <- format(uppervals, digits=digits, scientific=FALSE)
    rownames(confintMat) <- rownames(xresults)
    confintMat <- as.data.frame(confintMat)
    xresults <- cbind(xresults, confintMat)
  }
  
  ##print main table:
  cat("\n")
  cat("Estimation results:", "\n")
  cat("\n")
  print.data.frame(xresults, digits=digits)
  #format(xresults, digits=digits, scientific=FALSE)  
  #printCoefmat(xresults, digits=digits, signif.stars=FALSE, signif.legend=FALSE,
  #  P.values=FALSE)

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
##    digits        The minimum number of significant digits to be
##                  printed
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
gsummary.glm <- function(object, confint.level = 0.95, digits=6, ...)
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


  ##main table:
  ##-----------
  
  ##create:
  xnames <- c("coefs", "std.error", "t-value", "p-value")
  xresults <- matrix(NA, k, length(xnames))
  colnames(xresults) <- xnames
  rownames(xresults) <- names(coefs)
  xresults <- as.data.frame(xresults)
  
  ##fill:
  xresults[,"coefs"] <- format(coefs, digits=digits, scientific=FALSE)
  stderrors <- sqrt(diag(vcovmat))
  xresults[,"std.error"] <- format(stderrors, digits=digits, scientific=FALSE)
  tstats <- coefs/stderrors
  xresults[,"t-value"] <- format(tstats, digits=digits, scientific=FALSE)
  PvalTxt <- format(pt(abs(tstats), n-k, lower.tail=FALSE)*2, digits=digits,
    scientific=FALSE)
  xresults[,"p-value"] <- substr(PvalTxt, 1, digits+2)
    
  ##confidence intervals:
  if( !is.null(confint.level) ){
    confintMat <- matrix(NA, length(coefs), 2)
    confintMat <- as.data.frame(confintMat)
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
  print.data.frame(xresults, digits=digits)
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
