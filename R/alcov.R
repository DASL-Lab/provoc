#' Simplified Alcov Model Using Non-Negative Linear Regression
#'
#' A simplified Alcov model for linear regression without intercepts,
#' ensuring positivity of the coefficients. This function scales inputs and
#' applies a non-negative linear regression to estimate variant proportions.
#'
#' @param Y Vector of frequencies for each sample.
#' @param lmps Matrix or data frame of lineage definitions, similar to varmat.
#' @param muts Char vector of mutation names, used to select and order columns in lmps if it is a data frame.
#'
#' @return An "Alcov" object with coefficients from the linear regression model, representing the estimated proportions of variants. Object can be plotted.
#'
#' @examples
#' Y <- c(0.1, 0.2, 0.3, 0.4)
#' lmps <- matrix(c(0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1), nrow = 4, byrow = TRUE)
#' muts <- c("mut1", "mut2", "mut3")
#'
#' coef <- alcov(Y, lmps, muts)
#' print(coef)
#' 
#' If you want to plot the coefficients:
#' plot(coef, muts)
#'
#' @importFrom scales rescale
#' @importFrom nnls nnls
#' @export
alcov <- function(Y, lmps, muts) {
  # Ensure Y, lmps, and muts are properly aligned
  if (is.data.frame(lmps)) {
    lmps <- as.matrix(lmps[muts])
  }
  
  # Scale Y and lmps for better regression performance
  Y_scaled <- rescale(Y)
  lmps_scaled <- apply(lmps, 2, rescale)
  
  # Linear Regression without intercept using non-negative least squares
  model <- nnls(lmps_scaled, Y_scaled)
  
  # Extract and return coeffs
  alcov_coeffs <- coef(model)
  class(alcov_coeffs) <- "alcov"
  return(alcov_coeffs)
}


#' Plot Method for Alcov Coefficients
#'
#' @param coef Vector of coefficients returned by the alcov function.
#' @param muts Vector of mutation names, which must match the length of `coef`.
#' @importFrom ggplot2 ggplot geom_bar aes labs theme_minimal
#' @export
#' @method plot alcov
plot.alcov <- function(coef, muts) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 must be installed to use this function.")
  }
  
  df <- data.frame(Mutation = muts, Coefficient = coef)
  ggplot(df, aes(x = Mutation, y = Coefficient, fill = Mutation)) +
    geom_bar(stat = "identity") +
    labs(title = "Alcov Model Coefficients", x = "Mutation", y = "Coefficient") +
    theme_minimal()
}