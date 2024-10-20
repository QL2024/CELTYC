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
#'   A list object, containing 3 sublists of different DMCT sets.
#'   \itemize{
#'   \item dmct: A list object. Each entry is a vector of DMCTs altered in a cell type. The name of each entry is the cell type.
#'   \item spec: A list object. Each entry is a vector of DMCTs specific to a cell type. The name of each entry is the cell type.
#'   \item comm: A list object. Each entry is a vector of DMCTs common to 3 cell types. The name of each entry is the combination of 3 cell types.
#'   }
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
