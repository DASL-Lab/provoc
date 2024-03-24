
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
    b1, b2,
    newdata = NULL, type = NULL,
    dispersion = NULL, terms = NULL) {

    # TEMP example usage
    # -------
    # library(provoc)
    # data("Baaijens")
    
    # Baaijens$mutation <- parse_mutations(Baaijens$label)
    
    # res <- provoc(formula = cbind(count, coverage) ~ B.1.1.7 + B.1.617.2, data = Baaijens)
    
    # print(get_convergence(res))
    
    # b1 <- Baaijens[Baaijens$sra == unique(Baaijens$sra)[1], ]
    # b2 <- Baaijens[Baaijens$sra == unique(Baaijens$sra)[2], ]
    
    # predicted_values <- predict.provoc(res, b1, b2)
    # print(predicted_values)
    # -------
    # END OF TEMP

    if (!"provoc" %in% class(provoc_obj)) {
        stop("Object must be of class 'provoc'")
    }

    # b1_index <- match("b1", colnames(provoc_obj$variant))
    # b2_index <- match("b2", colnames(provoc_obj$variant))

    b1_index <- match(b1, provoc_obj$group)
    b2_index <- match(b2, provoc_obj$group)

    proportions <- as.numeric(provoc_obj$rho)
    variant_matrix <- attr(provoc_obj, "variant_matrix")
    #variant_matrix <- get_varmat(provoc_obj)

    if (any(!rownames(variant_matrix) %in% provoc_obj$variant)) {
        stop("Variant matrix does not match variants in results")
    }

    results <- proportions[b1_index] %*% variant_matrix[b2_index, ]
    #results <- proportions %*% variant_matrix[provoc_obj$variant, ]

    return(results)
}

resids.provoc <- function(provoc_obj, type = "deviance") {
    # TODO: Calculate deviance residuals
}
