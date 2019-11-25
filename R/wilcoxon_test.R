#' Wilcoxon Signed Rank and Rank-Sum Tests
#'
#' Performs a Wilcoxon Signed Rank test on a vector of data or two paired
#' vectors of data, or performs a Wilcoxon Rank-Sum test on two unpaired
#' vectors of data.
#' @usage wilcoxon_test(x, y = NULL, alternative = c("two.sided", "less", "greater"),
#' mu = 0, paired = FALSE, exact = NULL, correct = TRUE)
#'
#' @param x a numeric vector.
#' @param y an optional numeric vector.
#' @param alternative a string specifying the alternative hypothesis ("two.sided, "less", or "greater").
#' Only the first letter is necessary. Defaults to two sided.
#' @param mu a number representing a location or location difference to be tested for in the null hypothesis.
#' Defaults to 0.
#' @param paired a logical indicating whether to perform a paired or unpaired test.
#' Defaults to false.
#' @param exact a logical indicating whether to compute an exact p-value or an approximation.
#' @param correct A logical indicating whether to use a continuity correction or not.
#' @details If either y is not supplied or both x and y are supplied and paired is true,
#' then a Wilcoxon signed rank test is performed. The null hypothesis is that
#' the median of x (or x-y in the paired case) is equal to mu.
#'
#' If both x and y are supplied and paired is false, then a Wilcoxon rank-sum
#' (Mann-Whitney) test is performed. The null hypothesis is that the medians of
#' x and y differ by a location shift of mu. The alternative hypothesis "greater"
#' refers to x being shifted to the right of y.
#'
#' If exact is true, an exact p-value will be computed using the built-in
#' psignrank() or pwilcox() distribution functions. If exact is false, the
#' p-value is computed using a normal approximation. If exact is not supplied,
#' exact will default to true if the sample size is greater than 50 and there
#' are no ties. For the signed-rank test, there must also be no zeroes. Otherwise,
#' exact will be set to false.
#'
#' The correct parameter is only relevant if exact is false. If correct is true,
#' a continuity correction will be applied in order to perform the normal approximation.
#' The correction is a shifting of the test statistic by 0.5.
#'
#' Unlike the built in function wilcox.test(), this function does not support
#' calculating confidence intervals, or taking in a formula as an argument.
#'
#' @return A list of class "htest" consisting of the following:
#' \itemize{
#'   \item statistic - the value of the test statistic
#'   \item parameter - for this test, NULL
#'   \item p.value - the p-value
#'   \item null.value - the value of mu
#'   \item alternative - the alternative hypothesis
#'   \item method - the type of test performed
#'   \item data.name - the name(s) of the supplied data
#' }
#' @examples
#' wilcoxon_test(1:5)
#' @export


wilcoxon_test <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                          mu = 0, paired = FALSE, exact = NULL, correct = TRUE){
  alternative <- match.arg(alternative) # Input must correspond to one of the options
  data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  names(mu) <- "location shift"
  correction <- 0
  if(paired){
    if(is.null(y)){
      stop("'y' is missing for paired test")
    }
    else if(length(x) != length(y)){
      stop("'x' and 'y' must have the same length")
    }
    else{
      x <- x - y # Paired test is equivalent to one-sample test on x - y
      y <- NULL
    }
  }
  if(is.null(y)){ # If y is null or paired = T, perform a Wilcoxon signed rank test
    method <- "Wilcoxon signed rank test"
    if(!paired){
      data.name <- deparse(substitute(x))
      names(mu) <- "location"
    }
    x <- x - mu
    zero <- any(x == 0)
    if(zero){ # Remove zeroes, cannot perform exact test anymore
      x <- x[x != 0]
    }
    n <- length(x)
    ties <- n > length(unique(abs(x))) # If there are ties, cannot perform exact test anymore
    ranks <- rank(abs(x))
    signed_ranks <- ranks * sign(x)
    W <- sum(signed_ranks[signed_ranks > 0])
    names(W) <- "V"
    if(is.null(exact)){
      exact <- n < 50 & !ties & !zero # Use exact test for eligible small samples
    }
    if(exact & !ties & !zero){
      if(alternative == "two.sided"){
        if(W <= n * (n + 1) / 4){ # Check if W is smaller or greater than median
          p.value <- 2 * psignrank(W, n)
        }
        else{
          p.value <- 2 - 2 * psignrank(W - 1, n) # Subtract 1 because distribution is discrete
        }
      }
      else if(alternative == "less"){
        p.value <- psignrank(W, n)
      }
      else if(alternative == "greater"){
        p.value <- 1 - psignrank(W - 1, n)
      }
    }
    else if(exact){
      if(ties){
        warning("cannot compute exact p-value with ties")
      }
      if(zero){
        warning("cannot compute exact p-value with zeroes")
      }
      exact <- FALSE
    }
    if(!exact){ # If we can't use exact distribution, approximate with normal
      diff_W <- W - n * (n + 1) / 4
      var_W <- n * (n + 1) * (2 * n + 1) / 24
      numties <- table(ranks)
      var_W_reduction <- sum((numties ^ 3 - numties)) / 48 # Variance correction for ties
      var_W <- var_W - var_W_reduction
      if(correct){ # continuity correction
        if(alternative == "two.sided"){
          correction <- -sign(diff_W) / 2
        }
        else if(alternative == "less"){
          correction <- 0.5
        }
        else if(alternative == "greater"){
          correction <- -0.5
        }
        diff_W <- diff_W + correction
        method <- paste(method, "with continuity correction")
      }
      z <- diff_W / sqrt(var_W)
      if(alternative == "two.sided"){
        p.value <- 2 - 2 * pnorm(abs(z))
      }
      else if(alternative == "less"){
        p.value <- pnorm(z)
      }
      else if(alternative == "greater"){
        p.value <- 1 - pnorm(z)
      }
    }
  }
  else{ # If y is given and paired = F, perform a Wilcoxon rank-sum (Mann-Whitney) test.
    method <- "Wilcoxon rank sum test"
    x <- x - mu
    nx <- length(x)
    ny <- length(y)
    pooled <- c(x, y)
    ties <- nx + ny > length(unique(pooled))
    pooled_ranks <- rank(pooled)
    W <- sum(pooled_ranks[1:nx]) - nx * (nx + 1) / 2
    names(W) <- "W"
    if(is.null(exact)){
      exact <- nx < 50 & ny < 50 & !ties
    }
    if(exact & !ties){
      if(alternative == "two.sided"){
        if(W <= nx * ny / 2){ # Median is nx * ny / 2
          p.value <- 2 * pwilcox(W, nx, ny)
        }
        else{
          p.value <- 2 - 2 * pwilcox(W - 1, nx, ny)
        }
      }
      else if(alternative == "less"){
        p.value <- pwilcox(W, nx, ny)
      }
      else if(alternative == "greater"){
        p.value <- 1 - pwilcox(W - 1, nx, ny)
      }
    }
    else if(exact){
      warning("cannot compute exact p-value with ties")
      exact <- FALSE
    }
    if(!exact){
      diff_W <- W - nx * ny / 2
      var_W <- nx * ny * (nx + ny + 1) / 12
      numties <- table(pooled_ranks)
      var_W_reduction <- nx * ny * sum(numties ^ 3 - numties) / (12 * (nx + ny) * (nx + ny - 1))
      var_W <- var_W - var_W_reduction
      if(correct){ # continuity correction
        if(alternative == "two.sided"){
          correction <- -sign(diff_W) / 2
        }
        else if(alternative == "less"){
          correction <- 0.5
        }
        else if(alternative == "greater"){
          correction <- -0.5
        }
        diff_W <- diff_W + correction
        method <- paste(method, "with continuity correction")
      }
      z <- diff_W / sqrt(var_W)
      if(alternative == "two.sided"){
        p.value <- 2 - 2 * pnorm(abs(z))
      }
      else if(alternative == "less"){
        p.value <- pnorm(z)
      }
      else if(alternative == "greater"){
        p.value <- 1 - pnorm(z)
      }
    }
  }
  p.value <- unname(p.value)
  wilcoxon_list <- list(statistic = W, parameter = NULL, p.value = p.value,
                        null.value = mu, alternative = alternative,
                        method = method, data.name = data.name)
  class(wilcoxon_list) <- "htest"
  return(wilcoxon_list)
}
