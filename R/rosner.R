library(shades)
library(EnvStats)

#' Rosner's Test for Outliers
#'
#' Perform Rosner's generalized extreme Studentized deviate test for up to 
#' \eqn{k} potential outliers in a dataset, assuming the data without any outliers come 
#' from a normal (Gaussian) distribution. This is a modified version of a rosner
#' test from EnvStats package. 
#'
#' @param x numeric vector of observations. Missing (\code{NA}), undefined (\code{NaN}), 
#'   and infinite (\code{Inf}, \code{-Inf}) values are allowed but will be removed.  
#'   There must be at least 10 non-missing, finite observations in \code{x}.
#' @param k positive integer indicating the number of suspected outliers. The argument \code{k} 
#'   must be between 1 and \eqn{n-2} where \eqn{n} denotes the number of non-missing, finite 
#'   values in the argument \code{x}. The default value is \code{k=3}.
#' @param alpha numeric scalar between 0 and 1 indicating the Type I error associated with the 
#'   test of hypothesis. The default value is \code{alpha=0.05}.
#' @param warn logical scalar indicating whether to issue a warning (\code{warn=TRUE}; the default) 
#'   when the number of non-missing, finite values in \code{x} and the value of \code{k} are such 
#'   that the assumed Type I error level might not be maintained.
#' @param params string, indicating how to approach mean and sd estimation for the test
#' `iteratively` would be a default setting, like it's implemented in EnvStats
#'  where mean and as estimated each iteration, after removing the most extreme value
#' `extreme_params` would be replacing a mean and sd with a analogous location and
#'   scale parameters, after fitting an extreme distribution within a clone. Location and
#'   scale are estimated once, without updating it every iteration.
#' `remove_topx` would remove top x most extreme z-scores from the clone, leaving
#'  every potential outliers out, sinse we don't expect more than 2-3 specificities
#'  per clone.
#'  @param pval_dist what distribution type to use to calculate p value for R statistics
#'  `normal` pnorm function
#'  `t` pt function
#' @param remove_for_params for running rosner test in remove_topx setting in `params`
#' how much most extreme z-scores to remove, default is 3
#' 
#' @return A list with the results of the Rosner's test.
#'
#' @details This function performs Rosner's generalized extreme Studentized deviate test, which is 
#'   suitable for identifying multiple outliers in a dataset assumed to follow a normal distribution.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data <- rnorm(100)
#' rosnerTest(data, k = 3, alpha = 0.05)
#' }
#' 
#' @export
rosnerTest2 <- function(x, k = 3, alpha = 0.05, warn = TRUE, pval_dist = 't',
         params = 'iteratively', remove_for_params = 3){
  data.name <- deparse(substitute(x))
  if (!params %in% c('iteratively', 'remove_topx', 'extreme_params')){
    stop('params  must me either iteratively, remove_topx or extreme_params')
  }
  if (!is.numeric(x)) 
    stop("'x' must be a numeric vector")
  obs.num <- 1:length(x)
  if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
    x <- x[x.ok]
    obs.num <- obs.num[x.ok]
    warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
  }
  n <- length(x)
  if (n < 3) 
    stop("There must be at least 3 non-missing finite observations in 'x'")
  if (length(k) != 1 || !is.numeric(k) || !is.finite(k) || 
      k != round(k) || k < 1 || k > (n - 2)) 
    stop(paste("'k' must be a positive integer less than or equal to n-2,", 
               "where 'n' denotes the number of finite, non-missing observations in 'x'"))
  if (length(alpha) != 1 || !is.numeric(alpha) || !is.finite(alpha) || 
      any(alpha <= 0) || any(alpha >= 1)) 
    stop("'alpha' must be a numeric scalar greater than 0 and less than 1")
  if (warn) {
    if (k > 10 | k > floor(n/2)) {
      warning(paste("The true Type I error may be larger than assumed.", 
                    "Although the help file for 'rosnerTest' has a table with information", 
                    "on the estimated Type I error level,", "simulations were not run for k > 10 or k > floor(n/2).", 
                    sep = "\n"))
    }
    else {
      warn.conds <- (alpha > 0.01 & ((n >= 15 & n < 25 & 
                                        k > 2) | (n < 15 & k > 1))) | (alpha <= 0.01 & 
                                                                         (n < 15 & k > 1))
      if (warn.conds) 
        warning(paste("The true Type I error may be larger than assumed.", 
                      "See the help file for 'rosnerTest' for a table with information", 
                      "on the estimated Type I error level."))
    }
  }
  
  R <- rep(as.numeric(NA), k)
  x.vec <- rep(as.numeric(NA), k)
  obs.num.vec <- rep(as.numeric(NA), k)
  
  mean.vec <- rep(as.numeric(NA), k)
  sd.vec <- rep(as.numeric(NA), k)
  p_values <- rep(as.numeric(NA), k)  # Store p-values
  
  lambda <- rosnerTestLambda(n = n, k = 1:k, alpha = alpha)
  
  if (params == 'iteratively'){
    new.x <- x
    new.obs.num <- obs.num
    
    for (i in 1:k) {
      
      mean.vec[i] <- mean(new.x)
      sd.vec[i] <- sd(new.x)
      
      if (sd.vec[i] > 0) {
        abs.z = abs(new.x - mean.vec[i])/sd.vec[i]
        R[i] <- max(abs.z)
        index <- which(abs.z == R[i])[1]
        x.vec[i] <- new.x[index]
        obs.num.vec[i] <- new.obs.num[index]
        
        new.x <- new.x[-index]
        new.obs.num <- new.obs.num[-index]
      }
      else {
        R[i:k] <- NA
        break
      }
    }
  } else if (params == 'remove_topx'){
    
    if (length(x) <= remove_for_params+1){
      mean.x <- mean(x)
      sd.x <- sd(x)
    } else {
      subs <- remove_for_params-1
      
      mean.x <- mean(x[1:(length(x)-subs)]) 
      sd.x <- sd(x[1:(length(x)-subs)])
    }
    
    obs.num.vec <- order(x, decreasing = T)[1:k]
    x.vec <- x[obs.num.vec]
    abs.z = abs(x.vec - mean.x)/sd.x
    R <- abs.z[1:k]
    mean.vec <- rep(mean.x, k)
    sd.vec <- rep(sd.x, k)
    
    if (pval_dist == 'normal'){
      p_values <- pnorm(q = R / lambda, lower.tail = F)  
    } else if (pval_dist == 't') {
      df <- length(R) - (1:k) - 2
      df[df<=0] <- 1
      p_values <- pt(R, df = df, lower.tail = F) 
    }
  } 
  
  num.outlier.vec <- 1:k
  outlier <- R > lambda
  
  if (any(outlier)) {
    index <- max(num.outlier.vec[outlier], na.rm = TRUE)
    outlier[1:index] <- TRUE
  }
  
  out.df <- data.frame(num.outlier.vec - 1, mean.vec, sd.vec, 
                       x.vec, obs.num.vec, R, lambda, outlier, p_values)
  names(out.df) <- c("i", "Mean.i", "SD.i", "Value", "Obs.Num", 
                     "R.i+1", "lambda.i+1", "Outlier", "P.Values")
  distribution <- "Normal"
  dist.abb <- "norm"
  stat <- R
  names(stat) <- paste("R", 1:k, sep = ".")
  crit.value <- lambda
  names(crit.value) <- paste("lambda", 1:k, sep = ".")
  n.outliers <- sum(outlier, na.rm = TRUE)
  ret.list <- list(distribution = distribution, statistic = stat, 
                   sample.size = n, parameters = c(k = k), alpha = alpha, 
                   crit.value = crit.value, n.outliers = n.outliers, alternative = paste("Up to ", 
                                                                                         k, " observations are not\n", space(33), "from the same Distribution.", 
                                                                                         sep = ""), method = "Rosner's Test for Outliers", 
                   data = x, data.name = data.name, bad.obs = bad.obs, all.stats = out.df)
  oldClass(ret.list) <- "gofOutlier"
  
  ret.list
}


rosnerTestLambda <- function (n, k = 10, alpha = 0.05)
  {
    if (!is.numeric(n) || !all(is.finite(n)) || !all(n == round(n)) ||
        any(n < 3))
      stop("All values of 'n' must be positive integers greater than 2")
    if (!is.numeric(k) || !all(is.finite(k)) || !all(k == round(k)) ||
        any(k < 1) || any(k > (n - 2)))
      stop("All values of 'k' must be positive integers less than or equal to n-2")
    if (!is.numeric(alpha) || !all(is.finite(alpha)) || any(alpha <=
                                                            0) || any(alpha >= 1))
      stop("All values of 'alpha' must be greater than 0 and less than 1")
    arg.mat <- cbind.no.warn(n = as.vector(n), k = as.vector(k),
                             alpha = as.vector(alpha))
    for (i in c("n", "k", "alpha")) assign(i, arg.mat[, i])
    l <- k - 1
    p <- 1 - ((alpha/2)/(n - l))
    t.crit <- qt(p = p, df = n - l - 2)
    lambda <- t.crit * (n - l - 1)/sqrt((n - l - 2 + t.crit^2) *
                                          (n - l))
    names(lambda) <- NULL
    lambda
}

cbind.no.warn <- function (..., deparse.level = 1) 
  {
    oldopts <- options(warn = -1)
    on.exit(options(oldopts))
    base::cbind(..., deparse.level = deparse.level)
  }


#' Detect Outliers Using Extreme Value Theory
#'
#' This function applies an Extreme Value Theory (EVT) approach using the Gumbel distribution to detect outliers.
#' It estimates the distribution parameters and calculates p-values for each observation, based on different extreme
#' outlier detection strategies.
#'
#' @param x Numeric vector of values to test for extreme outliers.
#' @param type Character. Method to estimate extreme values. Options:
#'   - `'regular'` (default): Uses the estimated Gumbel distribution parameters.
#'   - `'reestimate_params'`: Iteratively re-estimates the distribution parameters, removing extreme values one at a time.
#'   - `'extreme_without_tail'`: Estimates distribution parameters after removing the three highest values.
#'
#' @return A numeric vector of p-values for each observation in `x`. Lower values indicate a higher likelihood of being an outlier.
#' If `reestimate_params` or `extreme_without_tail` is used, the function may return slightly different p-values based on the refined distribution.
#'
#' @examples
#' set.seed(123)
#' x <- c(rnorm(100, mean=50, sd=10), 120, 130, 140)  # Simulated data with outliers
#' pvals <- extreme_outlier_test(x, type="regular")
#' print(pvals)
#'
#' @export
extreme_outlier_test <- function(x, type='regular', double_loc_scale=F){
  
  extr <- extRemes::fevd(x = x, type = 'Gumbel')
  location <- extr$results$par["location"]
  scale <- extr$results$par["scale"]
  
  if (type == 'reestimate_params') { # might be an overkill
    pvalues_i <- extRemes::pevd(q = x, loc = location, scale = scale, type = 'Gumbel', lower.tail = F)
    
    k <- 5
    iter <- 0
    x_iter <- x
    while(pvalues_i[length(pvalues_i)] < extreme_alpha){
      iter <- iter + 1
      
      x_iter <- x_iter[-length(x_iter)]
      
      if (all(x_iter == 0)) {
        break
      }
      
      extr_i <- extRemes::fevd(x_iter, type = 'Gumbel')
      
      location_i <- extr_i$results$par["location"]
      scale_i <- extr_i$results$par["scale"]
      #shape_i <- extr_i$results$par["shape"]
      
      pvalues_i <- extRemes::pevd(q = x_iter, loc = location_i, scale = scale_i, type = 'Gumbel', lower.tail = F)
      
      if (iter >= k) {
        break
      }
    }
    pvalues <- extRemes::pevd(q = x, loc = location_i, scale = scale_i, type = 'Gumbel', lower.tail = F)
  } else if (type == 'extreme_without_tail') {
    
    x_new <- x[1:(length(x)-3)]
    
    extr_i <- extRemes::fevd(x_new, type = 'Gumbel')
    
    location_i <- extr_i$results$par["location"]
    scale_i <- extr_i$results$par["scale"]
    #shape_i <- extr_i$results$par["shape"]
    
    pvalues <- extRemes::pevd(q = x, loc = location_i, scale = scale_i, type = 'Gumbel', lower.tail = F)
  } else if (type == 'regular') {
    if (double_loc_scale){
      location <- location*2
      scale <- scale*2
    }
    pvalues <- extRemes::pevd(q = x, loc = location, scale = scale, type = 'Gumbel', lower.tail = F)
    #pvalues <- 1-extRemes::pevd(q = x, loc = location, scale = scale, shape = shape)
  } else {
    stop('extreme type can only by regular, reestimate_params or extreme_without_tail')
  }
  
  return(pvalues)
}
