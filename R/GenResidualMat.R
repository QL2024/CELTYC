#' Generate a residual matrix adjusted for cell type fractions from the original bulk DNAm data
#'
#' @description
#' This is a function to regress out estimated cell type fractions from the original DNAm data,
#' i.e. the beta matrix.
#'
#' @param beta.m
#' A beta matrix, with rows representing CpG sites and columns representing samples.
#'
#' @param standardize
#' A boolean value to determine whether to do standardization (zero-mean and unit-variance) for the residual matrix.
#' Default is TRUE.
#' @param estCTF.m
#' A matrix of cell type fractions, with columns representing cell types and rows representing samples.
#' @param ncores
#' The number of cores to use, i.e. at most how many child processes will be run simultaneously.
#'
#' @import parallel
#' @import stats
#'
#' @aliases GenResidualMat
#' @return
#' A (standardized if the parameter "standardized" is set to TRUE) residual matrix.
#'
#' @examples
#' \dontrun{
#' data("LIHC_data")
#' print(sd(LIHC_data$DNAm[1,]))
#' res.m <- GenResidualMat(LIHC_data$DNAm,estCTF.m = LIHC_data$estF,standardize = T,ncores = 40)
#' print(sd(res.m[1,]))
#' }
#'
#' @export
#'
#'
GenResidualMat <- function(beta.m,estCTF.m,standardize=T,ncores=4){
  res.o <- parallel::mclapply(1:nrow(beta.m),function(x){
    res.v <- lm(beta.m[x,]~estCTF.m)$res
    return(res.v)
  },mc.cores = ncores)
  res.m <- do.call(rbind,res.o)
  if(standardize==T){res.m <- t(scale(t(res.m)))}
  rownames(res.m) <- rownames(beta.m)
  colnames(res.m) <- colnames(beta.m)
  return(res.m)
}
