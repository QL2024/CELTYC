#' Extract DMCT sets from the output of CellDMC method
#'
#' @description
#' This function helps extract disease/phenotype associated differentially methylated CpG sites altered in one or more cell types (i.e. DMCTs)
#' from the output of the function  EpiDISH::CellDMC (\code{\link[EpiDISH]{CellDMC}}).
#'
#' It allows the users to retrieve DMCTs which are:
#' \enumerate{
#' \item Altered in each cell type (not necessarily the DMCTs specific to only one cell type).
#' \item Exclusively altered in each cell type.
#' \item Common to multiple cell types.
#'}
#'
#' @details
#' The retrieved DMCT sets can serve as the input of \code{\link[CELTYC]{DoCELTYC}} function.
#' Note that in order to do cell-type specific clustering properly, one should use a reasonable number of DMCTs.
#' Therefore, if there are not enough DMCTs for a particular cell type (fewer than the pre-set threshold parameter nDMCT),
#' getDMCT will not save that vector of DMCTs altered in that cell type.
#' Similarly, when users try to get DMCTs specific to a particular cell type (with the parameter "specific" set to TRUE),
#' if the number of the specific DMCTs is smaller than nSpDMCT,
#' getDMCT will not save that vector of DMCTs exclusively altered in that cell type.
#' When users try to get DMCTs common to nCT cell types (with the parameter "common" set to TRUE),
#' given a combination of nCT cell types,
#' if the number of the common DMCTs is smaller than nCmDMCT,
#' getDMCT will not save that vector of DMCTs common to these cell types.
#'
#' In general, these parameters (nDMCT, nSpDMCT, nCmDMCT) help users to keep only the DMCT sets containing a satisfying number of DMCTs.
#'
#'
#' @param celldmc.o A list object, which is the output of the function \code{\link[EpiDISH]{CellDMC}}
#' @param nDMCT
#' An integer value, which acts as a threshold value.
#' For a cell type, if the number of its DMCTs is smaller than nDMCT, the vector of its DMCTs will not be saved in the "dmct" entry of the returned list.
#'
#' @param specific Logical value, indicating whether to select DMCTs specific to each cell type or not. Default is TRUE.
#'
#' @param nSpDMCT
#' An integer value, which acts as a threshold value.
#' For a cell type, if the number of the DMCTs that are specific to that cell type is smaller than nSpDMCT, the vector of its specific DMCTs will not be saved in the "spec" entry of the returned list.
#'
#' @param common Logical value, indicating whether to select DMCTs which are common to different cell types. Default is TRUE.
#' @param nCT Integer value. The number of cell types for users to select common DMCTs.
#' @param nCmDMCT
#' An integer value, which acts as a threshold value.
#' For a combination of nCT cell types, if the number of the DMCTs common to these cell types is smaller than nCmDMCT, the vector of these common DMCTs will not be saved in the "comm" entry of the returned list.
#'
#' @return
#' A list with the following entries:
#' \itemize{
#' \item dmct: A list containing the vectors of DMCTs altered in different cell types.
#' Only DMCT vectors with a length >= nDMCT will be saved in this entry.
#' \item spec: A list containing the vectors of DMCTs exclusively altered in each of different cell types.
#' Only DMCT vectors with a length >= nSpDMCT will be saved in this entry. Note that this entry will be NULL if parameter "specific" is set to be FALSE.
#' \item comm: A list containing the vectors of DMCTs common to nCT cell types.
#' Only DMCT vectors with a length >= nCmDMCT will be saved in this entry. Note that this entry will be NULL if parameter "common" is set to be FALSE.
#' }
#'
#' @export
#'
GetDMCT <- function(celldmc.o,nDMCT=100,specific=T,nSpDMCT=100,common=T,nCT=3,nCmDMCT=100){
  ct.v <- colnames(celldmc.o[["dmct"]])[2:ncol(celldmc.o[["dmct"]])]
  dmct.lv <- list()
  for(ct in ct.v){
    dmct <- rownames(celldmc.o[["dmct"]])[which(celldmc.o[["dmct"]][,ct]!=0)]
    if(length(dmct)>=nDMCT){
      dmct.lv[[ct]] <-  rownames(celldmc.o[["dmct"]])[which(celldmc.o[["dmct"]][,ct]!=0)]}
  }
  ct.v <- names(dmct.lv) # select only cell types with more than nDMCT DMCTs
  message(paste0(sprintf("Identify more than %d DMCTs in these cell types: ",nDMCT),paste(ct.v,collapse = ", ")))
  if(length(dmct.lv)==0){
    stop(sprintf("None of the cell types have more than %d DMCTs.",nDMCT))
    }
  if(specific==T){ # get the DMCTs specific to each cell type
    spec.dmct.lv <- list()
      for(ct in ct.v){
        tmp.m <- celldmc.o[["dmct"]][dmct.lv[[ct]],setdiff(ct.v,ct)]
        num.v <- apply(tmp.m,1,function(x) length(which(x!=0)))
        other.dmct.v <- rownames(tmp.m)[which(num.v!=0)]
        cts.dmct.v <- setdiff(dmct.lv[[ct]],other.dmct.v)
        if(length(cts.dmct.v)>=nSpDMCT){
          spec.dmct.lv[[ct]] <- cts.dmct.v
        }else{
          message(sprintf("For %s, there are fewer than %d DMCTs that are exclusive to it.",ct,nSpDMCT))
        }
      }
  }else{
    spec.dmct.lv <- NULL
  }
  if(common==T){ # get the DMCTs shared by at least two cell types
      comm.dmct.lv <- list()
      if(nCT>length(ct.v)|nCT<2){
        if(length(ct.v)>2){
          stop(sprintf("Please input a value from the range of [2,%d] for the parameter nCM.",length(ct.v)))
        }else if(length(ct.v)==2){
          stop(sprintf("Please input 2 for the parameter nCM because there are only 2 cell types with more than %d DMCTs.",nDMCT))
        }else{
          stop(sprintf("It is not appropriate to select common DMCTs because there is only 1 cell type with more than %d DMCTs.",nDMCT))
        }
      }
      comb <- combn(names(dmct.lv),nCT,simplify = F)
      comm.dmct.lv <- lapply(comb, function(x) Reduce(intersect,dmct.lv[x],accumulate = F))
      names(comm.dmct.lv) <- lapply(comb,function(x) paste(x,collapse = "_"))
      rm.idx <- which(unlist(lapply(comm.dmct.lv,function(x) length(x)<nCmDMCT)))
      if(length(rm.idx)>0) comm.dmct.lv <- comm.dmct.lv[-rm.idx]
      if(length(comm.dmct.lv)==0){
        message(sprintf("For any combinations of the %d cell types, there are fewer than %d DMCTs that are common to them.",nCT,nCmDMCT))
      }
  }else{
    comm.dmct.lv <- NULL
  }
  return(list(dmct=dmct.lv,spec=spec.dmct.lv,comm=comm.dmct.lv))
}

