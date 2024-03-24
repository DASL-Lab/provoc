
#' Predict using Proportions of Variants of Concern
#'
#' Takes a named list with an estimate of the proportions and the
#' associated variant matrix, performs matrix multiplication to
#' predict outcomes, and returns results in the same order as the original data.
#'
#' @param provoc_obj Named list with `proportions` and `variant_matrix`.
#' @param newdata Not yet implemented.
#' @param type Not yet implemented.
#' @param se.fit Not yet implemented.
#' @param dispersion Not yet implemented.
#' @param terms Not yet implemented.
#' @param na.action Not yet implemented.
#' @return Predicted values in the same order as the input data.
#' @export
#' @examples
#' predicted_results <- predict(provoc_obj)
predict.provoc <- function(provoc_obj,
    newdata = NULL, type = NULL,
    dispersion = NULL, terms = NULL) {

    if (!"provoc" %in% class(provoc_obj)) {
        stop("Object must be of class 'provoc'")
    }

    proportions <- as.numeric(provoc_obj$rho)
    variant_matrix <- get_varmat(provoc_obj)
    if (any(!rownames(variant_matrix) %in% provoc_obj$variant)) {
        stop("Variant matrix does not match variants in results")
    }

    results <- proportions %*% variant_matrix[provoc_obj$variant, ]

    return(results)
}

resids.provoc <- function(provoc_obj, type = "deviance") {
    # TODO: Calculate deviance residuals
}


#' Autoplot for Proportions of Variants of Concern
#'
#' Generates a ggplot2 bar plot for the proportions of variants of concern from a `provoc` object.
#' Dynamically adjusts to handle both single and multiple samples within the same object.
#'
#' @param provoc_obj A `provoc` class object containing variant proportions and, optionally, multiple sample groups.
#' @import ggplot2
#' @return A ggplot object visualizing the variant proportions. If multiple samples are present,
#'         facets the plot for each sample.
#' @examples
#' autoplot.provoc(res)
#' @export
autoplot.provoc <- function(provoc_obj) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Please install.", call. = FALSE)
  }

  if (!"provoc" %in% class(provoc_obj)) {
    stop("Object must be of class 'provoc'")
  }
  
  # Prepare data for ggplot
  data_to_plot <- data.frame(variant = provoc_obj$variant,
                             proportion = provoc_obj$rho,
                             sample = ifelse(is.null(provoc_obj$group), "All", as.factor(provoc_obj$group)))
  
  p <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = variant, y = proportion, fill = variant)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Proportions of Variants of Concern",
                  x = "Variant", y = "Proportion") +
    ggplot2::scale_fill_viridis_d(begin = 0.3, end = 0.9, name = "Variant") # Using viridis for better color distinction
  
  # Facet the plot if there are multiple samples
  if (length(unique(data_to_plot$sample)) > 1) {
    p <- p + ggplot2::facet_wrap(~sample, scales = "free_x")
  }
  
  return(p)
}
