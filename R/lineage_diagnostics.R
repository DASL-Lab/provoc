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

#' Calculate the set differences for two character vectors
#' 
#' @param left,right Two character vectors
left_both_right <- function(left, right) {
    left_only <- setdiff(left, right)
    right_only <- setdiff(right, left)
    both <- intersect(left, right)

    list(left = left_only, mid = both, right = right_only)
}

#' Compare two lineage definition matrices
#' 
#' @param left_def,right_def Lineage definition matrice, such as those produced by \code{astronomize()} and \code{usher_barcodes()}.
#' @param col A named vector of colours, or a function with a single argument (n) that creates a colour scheme, such as \code{colorRampPalette()}. See Details and Examples.
#' @param  n_colours The number of colours to be used in the function that creates a diverging colour scheme.
#' @param ... Arguments passed on to \code{heatmap()}, such as \code{main} and \code{mar}.
#' 
#' @details A darker colour means that the mutation is present in that lineage, whereas a lighter colour means it's absent. 
#' 
#' Purple represents mutation/lineage combinations that are only present in the left definitions, and orange means it's only in the right.
#' 
#' The center square represents mutations and lineages that are present in both lineage definitions, with white indicating that the mutation is not present for that particular lineage and black indicating that it is present.
#' 
#' The colours can be a named vector, with names "NA_col", "left_0", "left_1", "NA_col", "right_0", "right_1", "neither", "left_1", "right_1", and "both". Alternatively, a diverging colour scale can be provided (the colour scale must hae a single argument that determines the number of colours).
#' 
#' @examples
#' 
#' lineages <- c("B.1.1.7", "P.1", "B.1.617.1", "B.1.617.2", "B.1.617.3")
#' left_def <- astronomize() |> filter_lineages(lineages[-1])
#' right_def <- usher_barcodes() |> filter_lineages(lineages[-2])
#' plot_lineage_defs2(left_def, right_def,
#'     main = "Constellations (without B.1.1.7) vs. Usher Barcodes (without P.1)")
#' 
#' my_cols <- colorRampPalette(c("dodgerblue", "white", "forestgreen"))
#' plot_lineage_defs2(left_def, right_def,
#'     main = "Constellations (without B.1.1.7) vs. Usher Barcodes (without P.1)",
#'     col = my_cols)
#' 
#' @seealso heatmap()
#' @export
plot_lineage_defs2 <- function(left_def, right_def,
    col = NULL, n_colours = 20, ...) {
    muts_lmr <- left_both_right(colnames(left_def), colnames(right_def))
    lins_lmr <- left_both_right(rownames(left_def), rownames(right_def))

    total_def <- matrix(ncol = length(unlist(muts_lmr)),
        nrow = length(unlist(lins_lmr)))
    colnames(total_def) <- unlist(muts_lmr)
    rownames(total_def) <- unlist(lins_lmr)
    defs <- c("left", "right", "neither", "left_only", "right_only", "both")
    # Left: 0 = NA, 1 = no, 2 = yes
    # Right: 3 = NA, 4 = no, 5 = yes
    # Both: 6 = neither, 7 = left_only,
    #     8 = right_only, 9 = both
    total_def <- matrix(0,
        ncol = length(unlist(muts_lmr)),
        nrow = length(unlist(lins_lmr)))
    colnames(total_def) <- unlist(muts_lmr)
    rownames(total_def) <- unlist(lins_lmr)
    # Lineages in left only
    total_def[lins_lmr$left, c(muts_lmr$mid, muts_lmr$left)] <-
        left_def[lins_lmr$left, c(muts_lmr$mid, muts_lmr$left)] + 1
    # Lineages in both mutations in left
    total_def[lins_lmr$mid, muts_lmr$left] <-
        left_def[lins_lmr$mid, muts_lmr$left] + 1
    # Lineages in right only
    total_def[lins_lmr$right, c(muts_lmr$mid, muts_lmr$right)] <-
        right_def[lins_lmr$right, c(muts_lmr$mid, muts_lmr$right)] + 4
    # Lineages in both, mutations in right
    total_def[lins_lmr$mid, muts_lmr$right] <-
        right_def[lins_lmr$mid, muts_lmr$right] + 4
    # Lineages in both, mutations in both
    total_def[lins_lmr$mid, muts_lmr$mid] <-
        left_def[lins_lmr$mid, muts_lmr$mid] +
        2 * right_def[lins_lmr$mid, muts_lmr$mid] +
        6

    if (is.null(col)) {
        col <- c(
            "NA_col" = "grey90",
            "left_0" = "darkorchid1",
            "left_1" = "darkorchid3",
            "right_0" = "chocolate1",
            "right_1" = "chocolate3",
            "neither" = "white",
            "both" = "grey40"
        )
    } else if (is.function(col)) {
        cols <- col(n_colours)
        col_breaks <- round(quantile(
            x = seq_len(n_colours),
            probs = c(0, 0.3, 0.5, 0.7, 1)
        ))
        col <- c(
            "NA_col" = "grey90",
            "left_0" = cols[col_breaks[2]],
            "left_1" = cols[col_breaks[1]],
            "right_0" = cols[col_breaks[4]],
            "right_1" = cols[col_breaks[5]],
            "neither" = cols[col_breaks[3]],
            "both" = "grey40"
        )
    }
    col_side_colours <- rep(col[c("left_1", "both", "right_1")],
        sapply(muts_lmr, length))
    row_side_colours <- rep(col[c("left_1", "both", "right_1")],
        sapply(lins_lmr, length))

    heatmap(total_def, Rowv = NA, Colv = NA,
        scale = "none",
        col = col[c("NA_col", "left_0", "left_1",
                "NA_col", "right_0", "right_1",
                "neither", "left_1", "right_1", "both")],
        ColSideColors = col_side_colours,
        RowSideColors = row_side_colours,
        revC = TRUE,
        ...)
}
