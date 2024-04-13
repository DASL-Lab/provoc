#' Plot a lineage definition matrix
#' 
#' @param lineage_defs A valid lineage definition matrix (rownames must be lineages)
#' @param col A function that takes argument \code{n} and returns a vector of colours of length \code{n}, such as \ode{colorRampPallette(c("white", "dodgerblue4"))}
#' @param main a main title for the plot.
#' 
#' @details Every other column in the display has a slight grey background to make it easier to read which mutation is in which lineage.
#' 
#' The final column summarises the number of mutations that are unique to a given variant, with the darker colour representing a variant with more unique mutations. The number in the brackets gives the exact number of unique mutations.
#' 
#' @examples
#' # Choose lineages that I know have a lot of shared mutations
#' lineage_defs <- astronomize() |>
#'     filter_lineages(c("B.1.1.7", "B.1.617.1", "B.1.617.2",
#'         "B.1.617.2+K417N", "B.1.617.3", "XBB.1"))
#' # XBB.1 has 1 shared mutation and 2 unique mutations
#' # None of B.1.1.7's mutations are shared
#' plot_lineage_defs(lineage_defs)
#' 
#' @export
plot_lineage_defs <- function(lineage_defs,
    col = colorRampPalette(colors = c("white", "dodgerblue4")),
    main = NULL) {

    ldef_processed <- lineage_defs |>
        filter_lineages(lineages = rownames(lineage_defs),
            shared_order = FALSE)

    uniques <- apply(
        X = ldef_processed,
        MARGIN = 2,
        FUN = function(x) sum(x) == 1
    )
    unique_counts <- apply(
        X = ldef_processed[, uniques],
        MARGIN = 1,
        FUN = sum
    )

    ldef_processed <- ldef_processed[, !uniques]

    # Turn every other column a light grey
    grey_cols <- 2 * (seq_len(ceiling(ncol(ldef_processed)/2))) - 1
    ldef_processed[, grey_cols] <- ldef_processed[, grey_cols] +
        max(unique_counts) + 3

    ldef_processed <- cbind(ldef_processed, "Unique" = unique_counts + 2)
    colnames(ldef_processed)[ncol(ldef_processed)] <- "Unique"

    rownames(ldef_processed) <- paste0(
        "(", unique_counts, ") ",
        rownames(ldef_processed)
    )

    first_cols <- c(0, 1)
    colours <- c(first_cols,
        col(max(ldef_processed) - length(first_cols) - 1),
        "grey95", 1)
    length(colours)
    max(ldef_processed)

    heatmap(ldef_processed, Rowv = NA, Colv = NA, scale = "none", 
        col = colours, main = main)
}
