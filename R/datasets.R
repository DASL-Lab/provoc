#' Facts about the lineages
#' 
#' Especially useful for checking the earliest/latest sequence date
#' 
#' @format A data frame
#' 
#' @source \url{https://nextstrain.org}
"lineage_facts"

#' Lists of mutations per lineage, as per Constellations
#' 
#' The names of mutations according to Pangolin's "Constellations" files. These are current as of October 2024. The constellation files for pre-Omicron lineages are not expected to change.
#' 
#' @format A data frame
#' 
#' @source \url{https://github.com/cov-lineages/constellations}
"constellation_lists"

#' Data from bioproject PRJNA741211 (Baaijens et al. (2022))
#' 
#' Includes 52 different samples in a single file, each sample being specified by the "sra" column. Data were processed according to a custom pipeline in https::/github.com/DASL-Lab/data-treatment-plant. The data were processed such that any mutation that had a frequency above 10% in any sample is ensured to be in *all* samples. That is, if C703T ever reached a frequency of 10% or more, then it is present in all samples (even if there were 0 observations) with the correct coverage at position 703.
#' 
#' The \code{label} column is in a proprietary format, with a "~" representing a SNP and "-" and "+" representing deletions and insertions, respectively. The code `parse_mutations(Baaijens$label)` will provide mutations in amino acid format, which is used throughout this package.
#' 
#' @references Baaijens, Jasmijn A., Alessandro Zulli, Isabel M. Ott, Ioanna Nika, Mart J. van der Lugt, Mary E. Petrone, Tara Alpert, Joseph R. Fauver, Chaney C. Kalinich, Chantal B. F. Vogels, Mallery I. Breban, Claire Duvallet, Kyle A. McElroy, Newsha Ghaeli, Maxim Imakaev, Malaika F. Mckenzie-Bennett, Keith Robison, Alex Plocik, Rebecca Schilling, Martha Pierson, Rebecca Littlefield, Michelle L. Spencer, Birgitte B. Simen, Ahmad Altajar, Anderson F. Brito, Anne E. Watkins, Anthony Muyombwe, Caleb Neal, Chen Liu, Christopher Castaldi, Claire Pearson, David R. Peaper, Eva Laszlo, Irina R. Tikhonova, Jafar Razeq, Jessica E. Rothman, Jianhui Wang, Kaya Bilguvar, Linda Niccolai, Madeline S. Wilson, Margaret L. Anderson, Marie L. Landry, Mark D. Adams, Pei Hui, Randy Downing, Rebecca Earnest, Shrikant Mane, Steven Murphy, William P. Hanage, Nathan D. Grubaugh, Jordan Peccia, Michael Baym, and Yale SARS-CoV-2 Genomic Surveillance Initiative. 2022. “Lineage Abundance Estimation for SARS-CoV-2 in Wastewater Using Transcriptome Quantification Techniques.” Genome Biology 23(1):236. doi: 10.1186/s13059-022-02805-9.
"Baaijens"
