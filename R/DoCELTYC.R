#' Perform cell-type-specific clustering
#'
#' @description
#' This function performs cell-type-specific clustering on input DNAm matrix.
#' @import ConsensusClusterPlus
#' @import r.jive
#'
#' @aliases DoCELTYC
#'
#' @param data.m A DNAm matrix on which to perform cell-type-specific clustering. A standardized residual matrix adjusted for cell type fractions is recommended here, but this is optional.
#' @param dmct.lv A list object where each entry is a character vector of DMCTs, i.e. phenotype associated differentially methylated CpGs altered in one or more cell types, which can be obtained from the output of function EpiDISH::CellDMC. 
#' The DMCT vector can be DMCTs for a specific cell type or a set of DMCTs shared among different cell types.
#' Please note that each entry should be assigned a name. The list of DMCTs can be retrieved by inputting the result of CellDMC into our function GetDMCT.
#'
#' @param method Character value. "jive": do JIVE analysis (with \code{\link[r.jive]{jive}}) on the input matrix first, and then do consensus clustering;
#' "consensus": do consensus clustering directly on the input matrix.
#' @param maxK Integer value. maximum cluster number to evaluate.
#' @param reps Integer value. number of subsamples.
#' @param pItem Numerical value. proportion of items to sample.
#' @param pFeature Numerical value. proportion of features to sample.
#' @param title Character value used as part of the output directory name for consensus clustering (for figures and tables).
#' @param clusterAlg 	Character value indicating cluster algorithm.
#' 'hc' for hierarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, or a function that returns a clustering.
#' @param distance 	Character value indicating the distance measure to be used. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
#' @param seed 	Optional numerical value. Sets random seed for reproducible results.
#' @param plot Character value. NULL - print to screen, 'pdf', 'png', 'pngBMP' for bitmap png, helpful for large datasets.
#' @param writeTable 	Logical value. TRUE - write ouput and log to csv.
#' @details
#' This function helps to do cell-type-specific clustering with 2 optional methods:
#' \itemize{
#' \item {jive}: By doing JIVE analysis first on the input matrices over different sets of DMCTs, one can obtain joint variation shared by the input DMCT sets
#'       and individual variation unique to each DMCT set.
#'       Note that when using "jive" method, there should be \strong{no overlapping} CpGs between every two entries of dmct.lv.
#'       Later consensus clustering will be performed on the joint variation matrix and different individual variation matrices derived from jive analysis.
#'       With this method the function returns a list with length equal to `length(dmct.lv) + 2`:
#'       The 1st entry contains the jive analysis result, and the rest contain the clustering results using joint variation (the 2nd entry) and individual variation matrices.
#' \item {consensus}: When using this method,
#' it is not necessary for the CpGs within each entry of "dmct.lv" to be completely non-overlapping. Consensus clustering will be performed directly on the input data matrix.
#' With this method the function returns a list with the same length as "dmct.lv", with each entry containing the consensus clustering results using different DMCT sets in "dmct.lv".
#' }
#'
#' @examples
#' \dontrun{
#'  data("LIHC_data")
#'  res.m <- GenResidualMat(LIHC_data$DNAm,estCTF.m = LIHC_data$estF,standardize = T,ncores = 40)
#'  test.results.l <- DoCELTYC(res.m,dmct.lv = LIHC_data$selDMCT[c("Lym")],method = "consensus",maxK = 3,title = "test")
#' }
#'
#' @return
#' A list of the same length as the input dmct.lv if the method is "consensus", or length(dmct.lv)+2 if method is "jive".
#' \itemize{
#' \item{For "consensus" method}: each entry contains the consensus clustering result using the DMCTs in each entry of "dmct.lv",
#' represented as a list of length maxK, and each internal list contains consensusMatrix (numerical matrix), consensusTree (hclust), consensusClass (consensus class assignments).
#' \item{For "jive" method}: the 1st entry is the jive analysis result, i.e. an object of class jive (\code{\link[r.jive]{jive}});
#' each subsequent entry contains the consensus clustering results, as in "consensus" method.
#' }
#'
#' @export
#'
#'
DoCELTYC <- function(data.m,dmct.lv,method="consensus",maxK=10,reps=1000,pItem=0.8,pFeature=1,
             title=NULL,clusterAlg="hc",distance="pearson",seed=123,
             plot="pdf",writeTable = T){
  common.dmct.lv <- list()
  results.l <- list()
  if(length(unique(names(dmct.lv)))<length(dmct.lv)) {
    stop("Error: each entry in dmct.lv needs a different name.")}
  for(i in 1:length(dmct.lv)){
    common.dmct.lv[[i]] <- intersect(rownames(data.m),dmct.lv[[i]])
  }
  names(common.dmct.lv) <- names(dmct.lv)

  if(method=="consensus"){
    message("Start doing consensus clustering on standardized residual matrix over input DMC list.")
    for(i in 1:length(common.dmct.lv)){
      results.l[[i]] = ConsensusClusterPlus::ConsensusClusterPlus(data.m[common.dmct.lv[[i]],],maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,
                                 title=paste0(title,"_",names(common.dmct.lv)[i]),clusterAlg=clusterAlg,
                                 distance=distance,seed=seed,plot=plot,writeTable = writeTable)
  # icl = ConsensusClusterPlus::calcICL(results,title=title,plot=pdf)
    }
    names(results.l) <- names(common.dmct.lv)
    return(results.l)
  }else if(method=="jive"){
    whole.lm <- list()
    for(i in 1:length(common.dmct.lv)){
      whole.lm[[i]] <- data.m[common.dmct.lv[[i]],]
    }
    names(whole.lm) <- names(common.dmct.lv)
    message("Start doing JIVE analysis.")
    jive.res <- r.jive::jive(whole.lm, showProgress = T)
    ## consensus clustering on the joint variation matrix derived from jive analysis
    tmp.m <- do.call(rbind,jive.res$joint)
    message("Start doing consensus clustering on the Joint variation matrix from jive analysis.")
    jive.JV.results = ConsensusClusterPlus::ConsensusClusterPlus(tmp.m,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,
                                              title=paste0("jive_JV_",title),clusterAlg=clusterAlg,
                                              distance=distance,seed=seed,plot=plot,writeTable = writeTable)
    ## consensus clustering on the individual variation matrices derived from jive analysis
    jive.IV.results.l <- list()
    message("Start doing consensus clustering on the Individual variation matrices from jive analysis.")
    for(i in 1:length(jive.res$individual)){
      jive.IV.results.l[[i]] = ConsensusClusterPlus::ConsensusClusterPlus(jive.res$individual[[i]],maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,
                                                    title=paste0("jive_IV_",title,"_",names(common.dmct.lv)[i]),clusterAlg=clusterAlg,
                                                    distance=distance,seed=seed,plot=plot,writeTable = writeTable)
    }
  jive.results.l <- list(jive_res=jive.res,JV=jive.JV.results)
  for(i in 1:length(jive.res$individual)){
    jive.results.l[[(i+2)]] <- jive.IV.results.l[[i]]
  }
  names(jive.results.l)[3:length(jive.results.l)] <- names(common.dmct.lv)
  return(jive.results.l)
  }else{stop("The method should be \"jive\" or \"consensus\".")}

}
