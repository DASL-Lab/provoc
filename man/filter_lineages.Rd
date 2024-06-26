% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/astronomize.R
\name{filter_lineages}
\alias{filter_lineages}
\title{Filter lineage_defs for specific lineages, keeping mutations that are present in at least one lineage}
\usage{
filter_lineages(
  lineage_defs = NULL,
  lineages = c("B.1.526", "B.1.1.7", "B.1.351", "B.1.617.2", "B.1.427", "B.1.429", "P.1"),
  return_df = FALSE,
  path = NULL,
  shared_order = TRUE
)
}
\arguments{
\item{lineage_defs}{The result of \code{astronomize()}. If NULL, tries to run \code{astronoimize}.}

\item{lineages}{Vector of lineage names (must be in \code{rownmaes(lineage_defs)}). Defaults to lineages circulating in 2021-2022.}

\item{return_df}{Should the function return a data frame? Note that returned df is transposed compared to lineage_defs. Default FALSE.}

\item{path}{Passed on to \code{astronomize} if \code{lineage_defs} is NULL.}

\item{shared_order}{Put shared mutations first? Default TRUE.}
}
\value{
A lineage definition matrix with fewer rows and columns than \code{lineage_defs}. If \code{return_df}, the columns represent lineage names and a \code{mutations} column is added.
}
\description{
Filter lineage_defs for specific lineages, keeping mutations that are present in at least one lineage
}
\details{
After removing some lineage, the remaining mutations might not be present in any of the remaining lineage. This function will remove mutations that no longer belong to any lineage.

\code{shared_order = TRUE} will result in the mutations that are present in the highest number of lineages to appear first. This is convenient for human inspection, but does not affect estimation.
}
\examples{
# After cloning the constellations repo
lineage_defs <- astronomize(path = "../constellations")
dim(lineage_defs)
lineage_defs <- filter_lineages(lineage_defs, c("B.1.1.7", "B.1.617.2"))
dim(lineage_defs) # rows and columns have changed

}
