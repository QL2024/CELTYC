#' LIHC DNAm beta matrix
#'
#'
#' A matrix, with rows for cancer-associated DMCs altered in Lymphocytes (identified with \code{\link[EpiDISH]{CellDMC}}),
#' columns for samples of cancer status
#'
#' @docType data
#'
#' @usage data(preCompData)
#'
#' @keywords datasets
#'
#' @name LIHC_DNAm
#'
"LIHC_DNAm"


#' Identified DMCTs for each cell type for LIHC DNAm data
#'
#'
#' A list, with each entry for DMCTs of each cell type
#'
#' @docType data
#'
#' @usage data(preCompData)
#'
#' @keywords datasets
#'
#' @name allDMCT
#'
#' @format
#' A list named "allDMCT" with 5 entries. Each entry is a vector of DMCTs for a cell type.
"allDMCT"


#' Estimated CTFs for LIHC DNAm data (for cancer status samples)
#'
#'
#' A matrix, with rows for samples and columns for cell types
#'
#' @docType data
#'
#' @usage data(preCompData)
#'
#' @keywords datasets
#'
#' @name LIHC_estF
#'
#' @format
#' A matrix
"LIHC_estF"

#' Phenotype information for LIHC DNAm data (for cancer status samples)
#'
#'
#' A matrix, with rows for samples and columns for different phenotypes
#'
#' @docType data
#'
#' @usage data(preCompData)
#'
#' @keywords datasets
#'
#' @name LIHC_pheno
#'
#' @format
#' A matrix
"LIHC_pheno"

#' Subsets of DMCTs for each cell type for LIHC DNAm data
#'
#'
#' A list, with each entry for DMCTs of each cell type.
#' For entries named "Lym", "Hep", "EC": these DMCTs are non-overlapped for these cell types.
#' For entry "common", these are DMCTs shared by 3 cell types.
#'
#' @docType data
#'
#' @usage data(preCompData)
#'
#' @keywords datasets
#'
#' @name selDMCT
#'
#'
"selDMCT"

#' A heatmap object of standardized residual matrix over lymphocyte DMCTs
#'
#'
#' A \code{\link[ComplexHeatmap]{Heatmap-class}} object, for a heatmap of
#' standardized residual matrix over lymphocyte specific DMCTs, annotated by CELTYC clustering results
#' using lymphocyte specific DMCTs, Hep specific DMCTs and EC specific DMCTs respectively.
#'
#' @docType data
#'
#' @usage data(preCompData)
#'
#' @keywords datasets
#'
#' @name clustHeatmap
#'
#' @format A \code{\link[ComplexHeatmap]{Heatmap-class}} object
"clustHeatmap"
