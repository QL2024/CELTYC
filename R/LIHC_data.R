#' All data related to the CELTYC analysis for LIHC samples in the tutorial.
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
#'   with rows labeling cancer-associated DMCTs in lymphocytes, hepatocytes and endothelial cells (identified with \code{\link[EpiDISH]{CellDMC}}),
#'   columns labeling cancer samples.}
#'   \item{pheno}{
#'   A data.frame, containing phenotype information for LIHC DNAm data,
#'   with rows labeling cancer samples and columns labeling different phenotypes.}
#'   \item{estF}{
#'   A matrix of estimated cell-type fractions for LIHC DNAm data, with rows labeling cancer samples and columns labeling cell types.}
#'   \item{allDMCT}{
#'   A list with 5 entries. Each entry is a vector of DMCTs for a cell type.
#'   }
#'   \item{selDMCT}{
#'   A list with each entry corresponding to a set of DMCTs.
#'   The entries labeled "Lym", "Hep", and "EC" contain DMCTs that
#'   are exclusive to individual cell types: lymphocytes, hepatocytes, and endothelial cells, respectively.
#'   The entry labeled "common" contain DMCTs that are shared among lymphocytes, hepatocytes, and endothelial cells.
#'   }
#'   \item{heatmap}{
#'   A heatmap object (\code{\link[ComplexHeatmap]{Heatmap-class}}),
#'   displaying the standardized residual matrix over lymphocyte-specific DMCTs,
#'   annotated by CELTYC clusters obtained using lymphocyte specific, hepatocyte specific, and endothelial cell (EC) specific DMCTs, respectively.
#'   }
#'   \item{jive_IV_lym}{
#'   A list containing the consensus clustering results (\code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}). It is obtained by running `DoCELTYC` using "jive" method.
#'   This list specifically stores the clustering results for the individual variation matrix for lymphocyte specific DMCTs.
#'   It has 3 entries. The k-th entry in it is the clustering results when cluster number is k: a list containing consensusMatrix (numerical matrix), consensusTree (hclust), consensusClass (consensus class assignments).
#'   }
#' }
#' @name LIHC_data
#'
"LIHC_data"
