% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuse.R
\name{fission}
\alias{fission}
\title{Un-fuse coco and lineage_defs.}
\usage{
fission(fused, sample = NULL)
}
\arguments{
\item{fused}{The result of \code{fuse(coco, lineage_defs)}}

\item{sample}{The name of the sample being used.}
}
\value{
A list containing coco and lineage_defs.
}
\description{
Fusion ensures that the mutation lists match and are in the correct order, but the two are separated for analysis.
}
\details{
This function is mainly used internally, but can be useful to check if which mutations are actually used in model fitting.
}
\keyword{internal}
