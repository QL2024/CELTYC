#' All data related to the CELTYC analysis for LIHC samples.
#'
#' This list object contains all data for LIHC samples that can be used in the vignette of this package.
#'
#' @docType data
#' @keywords LIHC data
#' @usage data("LIHC_data")
#'
#' @format A list with the following entries:
#' \describe{
#'   \item{DNAm}{
#'   A LIHC DNAm beta value matrix,
#'   with rows for cancer-associated DMCs altered in Lymphocytes (identified with \code{\link[EpiDISH]{CellDMC}}),
#'   columns for samples of cancer status.}
#'   \item{pheno}{
#'   A phenotype information matrix for LIHC DNAm data,
#'   with rows for samples of cancer status and columns for different phenotypes}
#'   \item{estF}{
#'   A matrix of estimated cell-type fractions for LIHC DNAm data, with rows for samples of cancer status and columns for cell types}
#'   \item{allDMCT}{
#'   A list with 5 entries. Each entry is a vector of DMCTs for a cell type.
#'   }
#'   \item{selDMCT}{
#'   A list with each entry for a set of DMCTs.
#'   For entries named "Lym", "Hep", "EC": These DMCTs are exclusive to each individual cell type.
#'   For entry "common", these are DMCTs shared by 3 cell types.
#'   }
#'   \item{heatmap}{
#'   A heatmap object (\code{\link[ComplexHeatmap]{Heatmap-class}}),
#'   for a heatmap of standardized residual matrix over lymphocyte specific DMCTs,
#'   annotated by CELTYC clustering results using lymphocyte specific DMCTs, Hep specific DMCTs and EC specific DMCTs respectively.
#'   }
#'   \item{jive_IV_lym}{
#'   A list for consensus clustering results (\code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}). It is obtained by running `DoCELTYC()` using "jive" method.
#'   This list specifically stores the results of doing clustering for individual variation matrix for lymphocyte specific DMCTs.
#'   This is a list of length 3 (which is the max cluster number). Each entry in it is a list containing consensusMatrix (numerical matrix), consensusTree (hclust), consensusClass (consensus class assignments).
#'   }
#' }
#' @name LIHC_data
#'
"LIHC_data"
