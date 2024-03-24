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
predict.provoc <- function(
    provoc_obj,
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


#' Plot Variant Proportions from a provoc Object
#'
#' Generates bar plots for the proportions of variants of concern (VoC) contained within a `provoc` object.
#' It can handle both single and multiple samples within the same `provoc` object.
#'
#' @param provoc_obj A `provoc` class object containing VoC proportion data and, optionally, multiple sample groups.
#'
#' @return Generates a bar plot or series of bar plots visualizing the proportion of each variant.
#'
#' @examples
#' plot.provoc(provoc_obj)
#'
#' @export
plot.provoc <- function(provoc_obj) {
    # Check if the object is of class 'provoc'
    if (!"provoc" %in% class(provoc_obj)) {
        stop("Object must be of class 'provoc'")
    }

    # Determine if there are multiple groups/samples
    multiple_samples <- !is.null(provoc_obj$group) && length(unique(provoc_obj$group)) > 1

    # Setting up the plot
    if (multiple_samples) {
        # Plot for multiple samples
        unique_groups <- unique(provoc_obj$group)
        par(mfrow = c(ceiling(length(unique_groups) / 2), 2)) # Adjusting the plotting area for multiple plots

        for (group in unique_groups) {
            group_data <- subset(provoc_obj, group == group)
            # Using rho as the height for the barplot
            barplot(
                height = group_data$rho,
                names.arg = group_data$variant,
                main = paste("Variant Proportions for Sample", group),
                xlab = "Variants",
                ylab = "Proportion",
                col = rainbow(nrow(group_data)), # Adding color
                las = 2
            ) # Rotating the x-axis labels for better readability
        }
    } else {
        # Plot for a single sample
        barplot(
            height = provoc_obj$rho,
            names.arg = provoc_obj$variant,
            main = "Variant Proportions",
            xlab = "Variants",
            ylab = "Proportion",
            col = rainbow(nrow(provoc_obj)), # Adding color
            las = 2
        ) # Rotating the x-axis labels for better readability
    }

    # Resetting the plot layout if it was changed
    if (multiple_samples) {
        par(mfrow = c(1, 1))
    }
}
