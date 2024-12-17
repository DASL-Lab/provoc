
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
