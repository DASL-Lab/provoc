#' Simplified Freyja Model Using Lasso Regression
#'
#' This function applies a Lasso regression to estimate proportions of variants
#' based on frequency and adjusting for sequence depth.
#'
#' @param mix Vector of frequencies (count divided by coverage) for each sample.
#' @param depths Vector of coverage for each sample.
#' @param df_barcodes Numeric matrix varmat, with rows as samples and columns as mutations.
#' @param muts Char vector of mutation names.
#' @param eps Very small number representing the strength for Lasso regression. Default is 1e-4.
#'
#' @return A "Freyja" object with coefficients from the Lasso regression, representing the estimated proportions of variants. Object can be plotted.
#'
#' @examples
#' mix <- c(0.1, 0.2, 0.3, 0.4)
#' depths <- c(10, 20, 30, 40)
#' df_barcodes <- matrix(c(0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1), nrow = 4, byrow = TRUE)
#' muts <- c("mut1", "mut2", "mut3")
#'
#' coef <- freyja(mix, depths, df_barcodes, muts)
#' print(coef)
#' 
#' If you want to plot the coefficients:
#' plot(coef, muts)
#'
#' @import glmnet
#' @export
freyja <- function(mix, depths, df_barcodes, muts, eps=1e-4) {
  # Adjust mutations based on coverage
  depth_adjustment <- log(depths + 1) / max(log(depths + 1))
  adjusted_mix <- mix * depth_adjustment
  
  # Adjusted to replicate depth adjustment for each mutation
  adjusted_barcodes <- t(t(df_barcodes) * depth_adjustment)

  # Prepare for glmnet
  x_matrix <- as.matrix(adjusted_barcodes)
  y_vector <- as.vector(adjusted_mix)
  
  # Initialize and fit the Lasso model
  lasso_model <- glmnet(x_matrix, y_vector, alpha = 1, lambda = eps, intercept = FALSE, lower.limits = 0)
  
  # Return the coeffs, excluding intercept
  # The intercept is the first element in the glmnet coefficient matrix, so we skip it
  freyja_coeffs <- coef(lasso_model)[-1]  # Removing intercept term which is included by default
  class(freyja_coeffs) <- "freyja"
  return(freyja_coeffs)
}


#' Plot Method for Freyja Coefficients
#'
#' @param coef Vector of coefficients returned by the freyja function.
#' @param muts Vector of mutation names, which must match the length of `coef`.
#' @importFrom ggplot2 ggplot geom_bar aes labs theme_minimal
#' @export
#' @method plot freyja
plot.freyja <- function(coef, muts) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 must be installed to use this function.")
  }
  
  df <- data.frame(Mutation = muts, Coefficient = coef)
  ggplot(df, aes(x = Mutation, y = Coefficient, fill = Mutation)) +
    geom_bar(stat = "identity") +
    labs(title = "Freyja Model Coefficients", x = "Mutation", y = "Coefficient") +
    theme_minimal()
}
