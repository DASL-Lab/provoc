% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lineage_diagnostics.R
\name{pairwise_lineage_plot}
\alias{pairwise_lineage_plot}
\title{Pairwise Jaccard similarity of all lineages from two separate definitions}
\usage{
pairwise_lineage_plot(left_def, right_def, prefix = c("", ""), ...)
}
\arguments{
\item{left_def, right_def}{Lineage definition matrices. Column names must be mutations in the same format, rownames can be arbitrary.}

\item{prefix}{Optional character vector of length 2. Add a prefix to the lineage names to better differentiate them. The first entry is taken to be the prefix for left_def, the second is for right_def.}

\item{...}{Options to be passed to \code{heatmap}, especially, main, col, margins}
}
\description{
Especially useful for arbitrary or made-up definitions, such as those from clustering of mutations.
}
