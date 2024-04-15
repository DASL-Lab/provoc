
#' Find position from amino acid string.
#'
#' @param aa A vector of mutations. aa:orf1a:I300V, ins:28215:3, del:27378:25, or C703T
#'
#' @details Finds the 0-indexed position of the amino acid. Intended to aid in adding coverage for mutations that weren't observed.
#' 
#' @examples
#' pos_from_aa("aa:orf1a:I300V") #1166
#' pos_from_aa("C703T") # 703
#'
#' @export
#' @keywords internal
pos_from_aa <- function(aa) {
    gcodes <- grepl("aa", aa)
    dels <- grepl("del", aa)
    inss <- grepl("ins", aa)
    points <- (!grepl(":", aa)) & grepl("[0-9][A-Z]$", aa)
    others <- !(gcodes | dels | inss | points)
    pos <- double(length(aa))

    orfs <- list( # note that these are 0-indexed
        "orf1a" = c(265, 13468),
        "orf1b" = c(13467, 21555),
        "S" = c(21562, 25384),
        "orf3a" = c(25392, 26220),
        "E" = c(26244, 26472),
        "M" = c(26522, 27191),
        "orf6" = c(27201, 27387),
        "orf7a" = c(27393, 27759),
        "orf7b" = c(2775, 27887),
        "orf8" = c(27893, 28259),
        "N" = c(28273, 29533),
        "orf10" = c(29557, 29674)
    )

    pos[gcodes] <- sapply(strsplit(aa[gcodes], ":"),
        function(x) {
            orfs[x[2]][[1]][1] + 3 * as.numeric(gsub("[^0-9]", "", x[3])) + 1
        })
    pos[dels] <- sapply(strsplit(aa[dels], ":"),
        function(x) {
            as.numeric(x[2])
        })
    pos[inss] <- as.numeric(sapply(strsplit(aa[inss], ":"),
        function(x) {
            as.numeric(x[2])
        }))
    pos[points] <- sapply(aa[points],
        function(x) {
            as.numeric(substr(x, 2, nchar(x) - 1))
        })
    pos[others] <- NA

    as.numeric(pos)
}

#' Get coverage at mutation positions.
#'
#' Currently requires "aa" notation.
#'
#' @param coverage A data frame with columns labelled "position" (0-indexed) and "coverage".
#' @param aa A vector of mutations. aa:orf1a:I300V, ins:28215:3, del:27378:25, or C703T
#' 
#' @details Finds the maximum coverage for the three positions of the amino acid. This ensures that the number of observed mutations can never be larger than the coverage.
#'
#' @export
#' @keywords internal
coverage_at_aa <- function(coverage, aa) {
    pos <- pos_from_aa(aa)
    get_three <- grepl("aa", aa)

    res <- sapply(which(!is.na(pos)), function(i) {
        ifelse(get_three[i],
            yes = max(coverage$coverage[coverage$position %in%
                        (pos[i]):(pos[i] + 2)]),
            no = coverage$coverage[coverage$position == pos[i]])
    })

    res[is.na(pos)] <- NA
    res
}

#' Add coverage of missing mutations to a data set
#'
#' An observation of 0 mutations out of 1 read is a different amount of information from 0 mutations out of 1,000 reads. This function adds the coverage for mutations used in \code{lineage_defs}, which benefits all models (but adds 0-inflation). 
#' 
#' @param coco The data set. 
#' @param coverage The coverage file, with position as the first column and coverage as the second.
#' @param mutation_list The mutations to be added (usually \code{colnames(lineage_defs)}).
#'
#' @return coco, but with added rows. The counts are assumed to be 0 for mutations that aren't in the data.
#' 
#' @details Currently assumes that the data has a column labelled "mutation" in amino acid notation.
#' 
#' @export
add_coverage <- function(coco, coverage, mutation_list) {
    new_muts <- mutation_list[!mutation_list %in% coco$mutation]
    cov_new_muts <- coverage_at_aa(coverage, new_muts)

    dplyr::bind_rows(coco,
        data.frame(mutation = new_muts, count = 0, coverage = cov_new_muts))
}
