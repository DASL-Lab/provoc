#' Create a lineage definition matrix from a list of lineages with their mutations
#'
#' The input should be a named list where the names represent lineages and the values are vectors of mutations. The output will be a valid lineage matrix.
#'
#' @param lineage_list A named list of vectors of mutations. Mutation names should be in the same format as those in your data, otherwise post-processing will be required.
#'
#' @return A lineage matrix (rownames are lineages, colnames are mutations, entry i,j is 1 if lineage i contains mutation j, 0 otherwise).
#' @export
lineage_defs_from_list <- function(lineage_list) {
    if (!is.list(lineage_list)) {
        stop("Input must be a list")
    } else if (is.null(names(lineage_list))) {
        stop("Input must be a named list.")
    } else if (any(sapply(lineage_list, function(x) class(x) != "character"))) {
        stop("Input must be a named list of character vectors")
    }

    all_mutations <- unique(unlist(lineage_list))
    lineage_defs <- as.data.frame(t(1 * sapply(lineage_list,
                function(x) all_mutations %in% x)))
    names(lineage_defs) <- all_mutations
    lineage_defs
}

#' Filter lineages active on a given date.
#'
#' Using NGSB data (see \code{?mutations_by_lineage}), checks that the earliest sequence in a lineage was observed by \code{start_date}. Optionally checks if the latest sequence in the lineage was observed after start_date. Optionally checks if the lineage was ever observed in Canada.
#'
#' @param lineage_names A character vector of lineage names. Must match the names in \code{mutations_by_lineage}.
#' @param start_date Earliest date in the study. Must be in ISO-8601 format (as should all dates with no exceptions).
#' @param check_after Check if there's an observed sequence after the start of your study? Set FALSE if the start date is recent.
#' @param check_canada Checks if the lineage was ever observed in Canada. Default FALSE.
#'
#' @details The code also ignores any + symbols and anything after them, so lineages such as B.1.617.2+K417N (Delta+) will be treated as B.1.617.2 (Delta).
#'
#' @examples
#' # BA.1 was only observed as of January 2021
#' extant_lineages(c("B.1.1.7", "B.1.617.2", "BA.1"), start_date = "2020-12-01")
#' 
#' # Subset lineage defs by date:
#' data(Baaijens)
#' b1 <- Baaijens[Baaijens$sra == Baaijens$sra[1], ]
#' max_date <- max(b1$date)
#' lineage_defs <- astronomize()
#' lins_to_check <- extant_lineages(rownames(lineage_defs), max_date)
#' lineage_defs <- filter_lineages(lineage_defs, lins_to check)
#' dim(lineage_defs)
#' 
#' @return A character vector.
extant_lineages <- function(lineage_names, start_date,
    check_after = TRUE, check_canada = FALSE) {

    start_date <- lubridate::ymd(start_date)
    lineage_names2 <- sapply(strsplit(lineage_names, "\\+"), `[`, 1)
    include_lineage <- logical(length(lineage_names2))

    for (i in seq_along(include_lineage)) {
        if (lineage_names2[i] %in% provoc::lineage_facts$pango_lineage) {
            lin_date <- lubridate::ymd(provoc::lineage_facts$earliest_date[
                lineage_facts$pango_lineage == lineage_names2[i]])
            include_lineage[i] <- lin_date <= start_date
        }

        if (check_after && include_lineage[i]) {
            late_date <- lubridate::ymd(provoc::lineage_facts$latest_date[
                lineage_facts$pango_lineage == lineage_names2[i]])
            include_lineage[i] <- include_lineage[i] & (late_date >= start_date)
        }
        if (check_canada & include_lineage[i]) {
            include_lineage[i] <- include_lineage[i] &
                !is.na(provoc::lineage_facts$Canada_count[
                    lineage_facts$pango_lineage == lineage_names2[i]])
        }
    }

    lineage_names[include_lineage]
}
