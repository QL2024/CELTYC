#' Perform cell-type-specific clustering
#'
#' @description
#' This function performs cell-type-specific clustering on input standardized residual matrix.
#' @import ConsensusClusterPlus
#' @import r.jive
#'
#' @param std.res.m A standardized residual matrix adjusted for cell type fractions.
#' @param method Character value. "jive": do JIVE on the input residual matrix first before doing clustering;
#' "consensus": do consensus clustering directly on the input residual matrix.
#' @param maxK integer value. maximum cluster number to evaluate.
#' @param reps integer value. number of subsamples.
#' @param pItem numerical value. proportion of items to sample.
#' @param pFeature numerical value. proportion of features to sample.
#' @param dmct.lv A list object. Each entry is a character vector of DMCs (CpG names) for a specific cell type
#' (or a set of DMCs shared by different cell types) to be used. Note that each entry should be given a name.
#' @param title character value for a part of the output directory name for consensus clustering (for figures and tables).
#' @param clusterAlg 	character value. cluster algorithm.
#' 'hc' hierarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, or a function that returns a clustering.
#' @param distance 	character value. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
#' @param seed 	optional numerical value. sets random seed for reproducible results.
#' @param plot character value. NULL - print to screen, 'pdf', 'png', 'pngBMP' for bitmap png, helpful for large datasets.
#' @param writeTable 	logical value. TRUE - write ouput and log to csv.
#' @details
#' This function helps to do cell type specific clustering with 2 optional methods: (1)"jive": By doing JIVE analysis first with input standardized residual matrices over different sets of DMCs unique to different cell types
#' or shared by multiple cell types, one can obtain joint variation shared by various cell types and individual variation unique to each cell type. Note that when using "jive" method, there should be **NO COMMON**
#' CpG probes in each entry of dmct.lv. Later consensus clustering will be performed on the joint variation matrix and different individual variation matricex derived from jive analysis. With this method the function returns
#' a list with length equal to `length(dmct.lv) + 2`: The 1st element is the jive analysis result, and the rest are the clustering results using joint variation and individual variation matrices. (2)"consensus": When using this method,
#' it is not necessary for the CpGs within each entry of dmct.lv to be completely non-overlapping.
#'
#' @return A list of length the same as length of dmct.lv or (length(dmct.lv)+2).
#' Each element is a list of length maxK, and each internal list contains consensusMatrix (numerical matrix), consensusTree (hclust), consensusClass (consensus class asssignments).
#' @export
#'
#'
DoCtsCluster <- function(std.res.m,method="consensus",maxK=10,reps=1000,pItem=0.8,pFeature=1,
             dmct.lv,title=NULL,clusterAlg="hc",distance="pearson",seed=123,
             plot="pdf",writeTable = T){
  common.dmct.lv <- list()
  results.l <- list()
  if(length(unique(names(dmct.lv)))<length(dmct.lv)) {
    stop("Error: each entry in dmct.lv needs a different name.")}
  for(i in 1:length(dmct.lv)){
    common.dmct.lv[[i]] <- intersect(rownames(std.res.m),dmct.lv[[i]])
  }
  names(common.dmct.lv) <- names(dmct.lv)

  if(method=="consensus"){
    message("Start doing consensus clustering on standardized residual matrix over input DMC list.")
    for(i in 1:length(common.dmct.lv)){
      results.l[[i]] = ConsensusClusterPlus::ConsensusClusterPlus(std.res.m[common.dmct.lv[[i]],],maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,
                                 title=paste0(title,"_",names(common.dmct.lv)[i]),clusterAlg=clusterAlg,
                                 distance=distance,seed=seed,plot=plot,writeTable = writeTable)
  # icl = ConsensusClusterPlus::calcICL(results,title=title,plot=pdf)
    }
    names(results.l) <- names(common.dmct.lv)
    return(results.l)
  }else if(method=="jive"){
    whole.lm <- list()
    for(i in 1:length(common.dmct.lv)){
      whole.lm[[i]] <- std.res.m[common.dmct.lv[[i]],]
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
