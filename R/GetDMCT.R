#' A function that returns a list of disease-associated DMCs altered in specific cell types (i.e. DMCTs)
#'
#' @description
#' This function employs a previously published algorithm CellDMC (https://doi.org/10.1038/s41592-018-0213-x)
#'
#'
#' @param beta.m
#' A beta value matrix with rows labeling the CpGs and columns labeling samples.
#'
#' @param pheno.v
#' A vector of phenotype. Both binary and continuous/oderinal phenotypes are allowed, but NA is not allowed in pheno.v.
#'
#' @param ctf.m
#' A matrix contains fractions of each cell-type. Each row labels a sample, with the same order of the columns in beta.m.
#' Each column labels a cell-type. Column names, which are the names of cell-types, are required. The rowSums of ctf.m should be 1 or close to 1.
#'
#' @param adjPMethod
#' The method used to adjust p values. The method can be any of method accepted by \code{\link{[stats]p.adjust}}.
#'
#' @param adjPThresh
#' A numeric value, default as 0.05. This is used to call DMCTs. For each cell-type respectively, the CpG with the adjusted p values less than this threshold will be reported as DMCTs (-1 or 1) in the 'dmct' matrix in the returned list.
#'
#' @param cov.mod
#' A design matrix from `model.matrix`, which contains other covariates to be adjusted. For example, input model.matrix(~ geneder, data = pheno.df) to adjust gender. Do not put cell-type fraction here!
#'
#' @param ncores
#' The number of cores to use, i.e. at most how many threads will run simultaneously.
#'
#' @import EpiDISH
#' @import stats
#'
#' @return
#' A list containing 2 lists. The 1st list is the returned results of function \code{\link[EpiDISH]{CellDMC}}.
#' The 2nd list is a list of identified DMCTs for each cell type.
#'
#'
#' @export
#'
GetDMCT <- function(beta.m,pheno.v,ctf.m,adjPMethod = "fdr",adjPThresh = 0.05,
                    cov.mod = NULL,ncores = 4){

  celldmc.o <- CellDMC(beta.m,pheno.v,ctf.m,mc.cores = ncores,cov.mod = cov.mod,
                       adjPMethod = adjPMethod,adjPThresh = adjPThresh)
  dmct.lv <- list() # DMCTs for each cell type
  for(ct in colnames(celldmc.o[["dmct"]])[2:ncol(celldmc.o[["dmct"]])]){
    dmct.lv[[ct]] <-  rownames(celldmc.o[["dmct"]])[which(celldmc.o[["dmct"]][,ct]!=0)]
  }
  return(list(celldmc_res=celldmc.o,dmct=dmct.lv))
}
