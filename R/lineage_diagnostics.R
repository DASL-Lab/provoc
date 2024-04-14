#' Plot a lineage definition matrix
#'
#' @param lineage_defs Lineage definition matrix, such as those produced by \code{astronomize()} and \code{usher_barcodes()}.Alternatively the output of \code{provoc()}, from which the lineage definition matrix will be extracted.
#' @param col A function that takes argument \code{n} and returns a vector of colours of length \code{n}, such as \code{colorRampPallette(c("white", "dodgerblue4"))}
#' @param main a main title for the plot.
#'
#' @details Every second column in the display has a slight grey background to make it easier to read which mutation is in which lineage.
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

    if (inherits(lineage_defs, "provoc")) {
        lineage_defs <- get_actual_defs(lineage_defs)
    }

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
#' @param left_def,right_def Lineage definition matrices, such as those produced by \code{astronomize()} and \code{usher_barcodes()}.Alternatively the output of \code{provoc()}, from which the lineage definition matrix will be extracted (one may be a lineage definition matrix and the other may be the result of fitting provoc, allowing for comparison of which lineages/mutations were actually used).
#' @param col A vector of three colours (left, middle, right). For a mutation present in \code{left_def} only, the colour will be the first colour in the vector, while a mutation not present will be 80% of the way between the left colour and the middle colour (similar for \code{right_def}). Optionally, the user can specify a vector of length 6, where the colours define left, middle, right, NA, neither, and both.
#' @param ... Arguments passed on to \code{heatmap()}, such as \code{main} for a title and \code{mar} to give the title space.
#'
#' @details A darker colour means that the mutation is present in that lineage, whereas a lighter colour means it's absent.
#'
#' The first colour in col represents mutation/lineage combinations that are only present in the left definitions, and the third color in col means it's only in the right.
#'
#' The center square represents mutations and lineages that are present in both lineage definitions, with white indicating that the mutation is not present for that particular lineage and black indicating that it is present.
#'
#' @examples
#'
#' lineages <- c("B.1.1.7", "P.1", "B.1.617.1", "B.1.617.2", "B.1.617.3")
#' left_def <- astronomize() |> filter_lineages(lineages[-1])
#' right_def <- usher_barcodes() |> filter_lineages(lineages[-2])
#' plot_lineage_defs2(left_def, right_def,
#'     main = "Constellations (without B.1.1.7) vs. Usher Barcodes (without P.1)")
#'
#' plot_lineage_defs2(left_def, right_def,
#'     main = "Constellations (without B.1.1.7) vs. Usher Barcodes (without P.1)",
#'     col = c(2, 0, 3))
#'
#' # See which mutations/lineages were actually used
#' lineage_defs <- astronomize() |>
#'     filter_lineages(c("B.1.617.1", "B.1.617.2", "B.1.617.2+K417N",
#'     "B.1.427", "B.1.429", "B.1.1.7"))
#' res <- provoc(count / coverage ~ ., data = b1,
#'     lineage_defs = lineage_defs)
#' # B.1.617.2 and B.1.617.2+K417N were "squashed" because they're identical.
#' # We can also see which mutations were *not* used in the analysis.
#' # They are the ones defined in the left_def but not present in right_def.
#' plot_lineage_defs2(lineage_defs, res)
#' 
#' @seealso heatmap()
#' @export
plot_lineage_defs2 <- function(left_def, right_def,
    col = c("#b2182b", "#f7f7f7", "#2166ac", "grey90", "white", "grey40"),
    ...) {

    if (inherits(left_def, "provoc")) {
        left_def <- get_actual_defs(left_def)
    }
    if (inherits(right_def, "provoc")) {
        right_def <- get_actual_defs(right_def)
    }

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

    if (any(col == 0)) {
        col[which(col == 0)] <- "white"
    }
    lmr <- colorRampPalette(col[1:3])(21)
    left_0 <- lmr[9]
    left_1 <- lmr[1]
    right_0 <- lmr[13]
    right_1 <- lmr[21]
    if (length(col) == 3) {
        na_col <- "grey90"
        neither <- "white"
        both <- "grey40"
    } else if (length(col) == 6) {
        na_col <- col[4]
        neither <- col[5]
        both <- col[6]
    } else {
        stop("col must be either 3 or 6 values")
    }

    col_side_colours <- rep(c(left_1, both, right_1),
        sapply(muts_lmr, length))
    row_side_colours <- rep(c(left_1, both, right_1),
        sapply(lins_lmr, length))

    heatmap(total_def, Rowv = NA, Colv = NA,
        scale = "none",
        col = c(na_col, left_0, left_1,
                na_col, right_0, right_1,
                neither, left_1, right_1, both),
        ColSideColors = col_side_colours,
        RowSideColors = row_side_colours,
        revC = TRUE,
        ...)
}