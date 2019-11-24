#' Wilcoxon Signed Rank and Rank-Sum Tests
#'
#' Performs a Wilcoxon Signed Rank test on a vector of data or two paired
#' vectors of data, or performs a Wilcoxon Rank-Sum test on two unpaired
#' vectors of data.
#'
#' @param x A number
#' @param y A number
#' @return The sum of \code{x} and \code{y}
#' @examples
#' add(1, 1)
#' add(10, 1)

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
    ties <- n > length(unique(x)) # If there are ties, cannot perform exact test anymore
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
  wilcoxon_list <- list(statistic = W, parameter = NULL, p.value = p.value,
                        null.value = mu, alternative = alternative,
                        method = method, data.name = data.name)
  class(wilcoxon_list) <- "htest"
  return(wilcoxon_list)
}
