# TODO: Better organization of functions. The plotting functions and the similarity calculations should be separate.

#' Plot a lineage definition matrix
#'
#' @param lineage_defs Lineage definition matrix, such as those produced by \code{astronomize()} and \code{usher_barcodes()}. Alternatively the output of \code{provoc()}, from which the lineage definition matrix will be extracted.
#' @param col A function that takes argument \code{n} and returns a vector of colours of length \code{n}, such as \code{colorRampPallette(c("white", "dodgerblue4"))}
#' @param ... Further arguments to be passed to \code{heatmap}, especially \code{main} and \code{margins}.
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
    col = colorRampPalette(colors = c("white", "dodgerblue4")), main = "Lineage Definitions", ...) {

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
        col = colours, main = main, ...)
}

#' Calculate the set differences for two character vectors
#'
#' @param left,right Two character vectors
#' @keywords internal
left_both_right <- function(left, right) {
    left_only <- setdiff(left, right)
    right_only <- setdiff(right, left)
    both <- intersect(left, right)

    list(left = left_only, mid = both, right = right_only)
}

#' Compare two lineage definition matrices
#' 
#' TODO: Check where this is used, make sure it's as useful as possible, give a more informative name.
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
        if (!inherits(left_def, "matrix")) {
            left_def <- total_lineage_defs(left_def, summarise = "or")
        }
    }
    if (inherits(right_def, "provoc")) {
        right_def <- get_actual_defs(right_def)
        if (!inherits(right_def, "matrix")) {
            right_def <- total_lineage_defs(right_def, summarise = "or")
        }
    }
    left_def <- as.matrix(left_def)
    right_def <- as.matrix(right_def)


    mutsets <- left_both_right(colnames(left_def), colnames(right_def))
    linsets <- left_both_right(rownames(left_def), rownames(right_def))

    total_def <- matrix(ncol = length(unlist(mutsets)),
        nrow = length(unlist(linsets)))
    colnames(total_def) <- unlist(mutsets)
    rownames(total_def) <- unlist(linsets)
    defs <- c("left", "right", "neither", "left_only", "right_only", "both")
    # Left: 0 = NA, 1 = no, 2 = yes
    # Right: 3 = NA, 4 = no, 5 = yes
    # Both: 6 = neither, 7 = left_only,
    #     8 = right_only, 9 = both
    total_def <- matrix(0,
        ncol = length(unlist(mutsets)),
        nrow = length(unlist(linsets)))
    colnames(total_def) <- unlist(mutsets)
    rownames(total_def) <- unlist(linsets)

    # Lineages in left only
    if (length(linsets$left) > 0 && length(c(mutsets$mid, mutsets$left)) > 0)
        total_def[linsets$left, c(mutsets$mid, mutsets$left)] <-
            left_def[linsets$left, c(mutsets$mid, mutsets$left)] +
            1

    # Lineages in both mutations in left
    if (length(linsets$mid) > 0 && length(mutsets$left) > 0)
        total_def[linsets$mid, mutsets$left] <-
            left_def[linsets$mid, mutsets$left] +
            1

    # Lineages in right only
    if (length(linsets$right) > 0 && length(c(mutsets$mid, mutsets$right)) > 0)
        total_def[linsets$right, c(mutsets$mid, mutsets$right)] <-
            right_def[linsets$right, c(mutsets$mid, mutsets$right)] +
            4

    # Lineages in both, mutations in right
    if (length(linsets$mid) > 0 && length(mutsets$right) > 0)
        total_def[linsets$mid, mutsets$right] <-
            right_def[linsets$mid, mutsets$right] +
            4

    # Lineages in both, mutations in both
    if (length(linsets$mid) > 0 && length(mutsets$mid) > 0)
        total_def[linsets$mid, mutsets$mid] <-
            left_def[linsets$mid, mutsets$mid] +
            2 * right_def[linsets$mid, mutsets$mid] +
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
        sapply(mutsets, length))
    row_side_colours <- rep(c(left_1, both, right_1),
        sapply(linsets, length))

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

#' Adds a list of lineage definition matrices together.
#' 
#' Entries represent the total number of times that lineage/mutation combination were used during fitting. Useful for checking whether a mutation was not present in a given sample.
#' 
#' @param ldef_list A list of lineage definition matrices, such as those extracted by \code{get_actual_defs}.
#' @keywords internal
total_lineage_defs <- function(ldef_list, summarise = "add") {
    all_muts <- sapply(ldef_list, colnames) |> as.character() |> unique()
    all_lins <- sapply(ldef_list, rownames) |> as.character() |> unique()

    total_def <- matrix(0, ncol = length(all_muts),
        nrow = length(all_lins))
    colnames(total_def) <- all_muts
    rownames(total_def) <- all_lins

    for (i in seq_along(ldef_list)) {
        these_muts <- colnames(ldef_list[[i]])
        these_lins <- rownames(ldef_list[[i]])
        if (summarise == "add") {
            total_def[these_lins, these_muts] <-
                total_def[these_lins, these_muts] +
                ldef_list[[i]][these_lins, these_muts]
        } else {
            total_def[these_lins, these_muts] <-
                total_def[these_lins, these_muts] |
                ldef_list[[i]][these_lins, these_muts]
        }
    }

    zero_cols <- which(apply(total_def, 2, sum) == 0)
    length(zero_cols)

    total_def
}


#' Summarise coverage by lineage/mutation combo
#' 
#' Useful for finding whether some calls were based on less information
#' 
#' @param provoc_obj Result of `provoc()`.
#' @param fun The function used to summarise the coverage information.
#' @param ... Further arguments passed onto \code{fun()}, such as \code{prob = 0.5} to find the median when fun = \code{quantile}.
#' 
#' @return A matrix with the total column and row names of all lineage defs, and entries indicating the number of analyses which used that mutation in that lineage.
coverage_by_lineage_defs <- function(provoc_obj, fun = mean, ...) {
    fused <- attributes(provoc_obj)$internal_data

    all_muts <- unique(fused$mutation)
    all_lins <- colnames(fused)[startsWith(colnames(fused), "lin")]

    total_def <- matrix(0, ncol = length(all_muts),
        nrow = length(all_lins))
    colnames(total_def) <- all_muts
    rownames(total_def) <- all_lins

    for (mut in all_muts) {
        for (lin in all_lins) {
            cond <- (fused$mutation == mut) &
                (fused[, lin] == 1)
            if (sum(cond) == 0) {
                total_def[lin, mut] <- NA
            } else {
                total_def[lin, mut] <- fun(
                    fused$coverage[fused$mutation == mut &
                            fused[, lin] == 1],
                    ...
                )
            }
        }
    }

    rownames(total_def) <- gsub("lin_", "", rownames(total_def))
    total_def
}

#' Visualize how lineage definitions were used in analysis
#' 
#' TODO: Is this a duplicate of plot_lineage_defs() above??? Check with real data and choose the one I like best.
#' 
#' @param provoc_obj The result of calling \code{provoc}.
#' @param type "used" for the total number of times a mutation was used in the analysis and "coverage" for information about the coverage. Default \code{"coverage"}.
#' @param fun The function used to summarise coverage (when type = "coverage"). Default \code{mean}.
#' @param ... Further arguments passed on to \code{fun}
#' @param col Colours to be used in plotting. Default \code{hcl.colors(n = 21, palette = "Dark Mint", rev = TRUE)}. Single-hue sequential colour palettes recommended. 
#' @param main A main title for the plot.
#' @param margins The margins of the plot, as used by \code{heatmap()}.
#' 
#' @examples
#' 
#' data(Baaijens)
#' b2 <- Baaijens [Baaijens$sra %in% unique(Baaijens$sra)[1:30], ]
#' b2$mutations <- parse_mutations(b2$label)
#' 
#' lineage_defs <- astronomize() |>
#'     filter_lineages(c("B.1.617.1", "B.1.617.2", "B.1.617.2+K417N",
#'         "B.1.427", "B.1.429", "B.1.1.7"))
#' res <- provoc(
#'     formula = count / coverage ~ .,
#'     lineage_defs = lineage_defs,
#'     data = b2, 
#'     by = "sra",
#'     bootstrap_samples = 0)
#' plot_actual_defs(res, type = "coverage", fun = max, na.rm = TRUE)
#' @export
plot_actual_defs <- function(provoc_obj,
    type = "coverage", fun = mean, ..., main = NULL,
    col = hcl.colors(n = 21, palette = "Dark Mint", rev = TRUE),
    margins = c(12, 10)) {
    if (type == "used") {
        total_def <- total_lineage_defs(get_actual_defs(provoc_obj))
    } else if (type == "coverage") {
        total_def <- coverage_by_lineage_defs(provoc_obj, fun = fun, ...)
    }

    rowmax <- apply(total_def, 1, max, na.rm = TRUE) |> round(2)
    rowmed <- apply(
        X = total_def,
        MARGIN = 1,
        FUN = function(x) quantile(x[x > 0], 0.5, na.rm = TRUE)
    ) |>
        round(2)
    rowmin <- apply(total_def, 1, function(x) min(x[x > 0], na.rm = TRUE)) |>
        round(2)

    rownames(total_def)[1] <- paste0(rownames(total_def)[1],
        "\n(min, median, max) = ")
    rownames(total_def) <- paste0(
        rownames(total_def),
        "\n(", rowmin, ", ", rowmed, ", ", rowmax, ")"
    )

    if (is.null(main)) {
        if (type == "coverage") {
            main <- paste0(as.character(substitute(fun)),
                " of coverage within each sample.")
        } else {
            main <- "Total times each definition was used."
        }
    }

    heatmap(total_def, Rowv = NA, Colv = NA, scale = "none",
        revC = TRUE, col = col, main = main, margins = margins)
}

#' Pairwise Jaccard similarity of all lineages from two separate definitions
#' 
#' Especially useful for arbitrary or made-up definitions, such as those from clustering of mutations.
#' 
#' TODO: Is this a duplicate of plot_lineage_defs2() above???
#' 
#' @param left_def,right_def Lineage definition matrices. Column names must be mutations in the same format, rownames can be arbitrary.
#' @param prefix Optional character vector of length 2. Add a prefix to the lineage names to better differentiate them. The first entry is taken to be the prefix for left_def, the second is for right_def.
#' @param ... Options to be passed to \code{heatmap}, especially, main, col, margins
#' 
#' @export
pairwise_lineage_plot <- function(left_def, right_def,
    prefix = c("", ""), ...) {
    shared_muts <- intersect(colnames(left_def), colnames(right_def))
    # Rows are left_def, cols are right_def
    lineage_similarity <- matrix(0,
        nrow = nrow(left_def),
        ncol = nrow(right_def))
    rownames(lineage_similarity) <- paste0(prefix[1],
        rownames(left_def))
    colnames(lineage_similarity) <- paste0(prefix[2],
        rownames(right_def))
    for (lin1 in seq_along(rownames(lineage_similarity))) {
        for (lin2 in seq_along(colnames(lineage_similarity))) {
            l2 <- right_def[rownames(right_def)[lin2], shared_muts]
            l1 <- left_def[rownames(left_def)[lin1], shared_muts]
            intersection <- as.logical(l1) & as.logical(l2)
            lineage_similarity[lin1, lin2] <- sum(intersection) /
                length(shared_muts)
        }
    }
    heatmap(lineage_similarity,
        Rowv = NA, Colv = NA, scale = "none", ...)
}
