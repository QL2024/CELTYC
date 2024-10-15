## ----chunk1-------------------------------------------------------------------
knitr::opts_chunk$set(crop = NULL)
library(CELTYC)
data("LIHC_data") # load the data used in the following analysis

## ----chunk2, eval=F, echo=T---------------------------------------------------
#  library(EpiSCORE)
#  avDNAm.m <- constAvBetaTSS(bmiq.m,type="450k") #computing the average DNAm over a window 200bp upstream of TSS, or if not available over 1st exon probes.
#  estF.o <- wRPC(avDNAm.m,ref=mrefLiver.m,useW=TRUE,wth=0.4,maxit=500) # Estimating cell-type fraction with a DNAm reference matrix for liver via weighted robust linear regression
#  estF.m <- estF.o$estF

## ----chunk3, eval=F, echo=T---------------------------------------------------
#  phe.v <- pheno.df$cancer # the cancer status for samples
#  phe.v[which(phe.v=="Cancer")] <- 1
#  phe.v[which(phe.v=="Normal")] <- 0
#  phe.v <- as.numeric(phe.v)
#  sex.v <- pheno.df$gender
#  sex.v[which(sex.v=="MALE")] <- 1
#  sex.v[which(sex.v=="FEMALE")] <- 2
#  sex.v <- as.numeric(sex.v)
#  age.v <- pheno.df$age_at_initial_pathologic_diagnosis
#  na.idx <- which(is.na(age.v)|is.na(sex.v)) # the running of CellDMC does not allow NA in input
#  
#  cov.mod.use <- model.matrix(~sex.v[-na.idx]+age.v[-na.idx]) # exclude the confounding effect of age and sex
#  celldmc.o <- CellDMC(bmiq.m[,-na.idx],pheno.v=phe.v[-na.idx],ctf.m = estF.m[-na.idx,],
#                       mc.cores = 40,ctf.m = estF.m[-na.idx,],cov.mod = cov.mod.use)
#  dmct.lv <- list() # save DMCTs for each cell type
#  for(ct in colnames(celldmc.o[["dmct"]])[2:ncol(celldmc.o[["dmct"]])]){
#      dmct.lv[[ct]] <-  rownames(celldmc.o[["dmct"]])[which(celldmc.o[["dmct"]][,ct]!=0)]
#  }
#  dmct.res <- list(celldmc_res=celldmc.o,dmct=dmct.lv)
#  

## ----chunk4, message=F--------------------------------------------------------
library(SuperExactTest)
supertest.o <- supertest(LIHC_data$allDMCT,n=403197)
plot(supertest.o,"landscape",sort.by = "size",color.scale.cex = 0.7,overlap.size.cex = 0.6,cex=0.7)

## ----chunk5,echo=T------------------------------------------------------------
res.m <- GenResidualMat(LIHC_data$DNAm,estCTF.m = LIHC_data$estF,standardize = T,ncores = 40) # regress out CTFs from DNAm matrix and get standardized residual matrix
CELTYC.results.l <- DoCELTYC(res.m,method = "consensus",maxK = 3,dmct.lv = LIHC_data$selDMCT[c("Lym","Hep","EC")],title = "cluster-test") # do consensus clustering on standardized residual matrix over DMCTs specific to lymphocytes, hepatocytes and endothelial cells

## ----chunk6,message=F---------------------------------------------------------
library(ComplexHeatmap)
lym.clust.v <- CELTYC.results.l$Lym[[3]]$consensusClass
hep.clust.v <- CELTYC.results.l$Hep[[3]]$consensusClass
ec.clust.v <- CELTYC.results.l$EC[[2]]$consensusClass # if K=3, the 3rd cluster contains few samples
print("Compare the lym-clusters and hep-clusters:")
table(lym.clust.v,hep.clust.v)
print("Compare the lym-clusters and EC-clusters:")
table(lym.clust.v,ec.clust.v)
print("Compare the hep-clusters and EC-clusters:")
table(hep.clust.v,ec.clust.v)

ComplexHeatmap::draw(LIHC_data$heatmap)

## ----chunk7-------------------------------------------------------------------
library(survival)

# for clustering results using DMCTs for different cell types:
clust.l <- list(Lym=lym.clust.v,Hep=hep.clust.v,EC=ec.clust.v)
surv.res.l <- list()
for(m in 1:3){
  clust.v <- clust.l[[m]]
  chi.surv.p.v <- vector()
  tmp.v <- vector()
  count <- length(unique(na.omit(clust.v)))
  cl <- sort(unique(clust.v))

  # generate pairwise P values comparing survival probability between different clusters
  for(i in 1:(count-1)){
    for(j in (i+1):count){
      tmp.time <- LIHC_data$pheno$OS.time[clust.v %in% c(cl[i],cl[j])]
      tmp.event <- LIHC_data$pheno$event[clust.v %in% c(cl[i],cl[j])]
      tmp.clust.v <- clust.v[clust.v %in% c(cl[i],cl[j])]
      surv.o <- Surv(time = tmp.time,event = tmp.event)
      cox.o <- coxph(surv.o~tmp.clust.v)
      chisq.pval <- pchisq(cox.o$score,1,lower.tail = F)
      chi.surv.p.v <- c(chi.surv.p.v,chisq.pval)
      tmp.v <- c(tmp.v,paste0("cl ",cl[i],",",cl[j]))
    }
  }  
  names(chi.surv.p.v) <- tmp.v
  # generate overall KM-curve plot to see the difference of survival rate between clusters
  surv.o <- Surv(time = LIHC_data$pheno$OS.time,event = LIHC_data$pheno$event)
  survfit.o <- survfit(surv.o ~ clust.v) 
  surv.res.l[[m]] <- list()
  surv.res.l[[m]][["pair-Pval"]] <- chi.surv.p.v
  surv.res.l[[m]][["survfit"]] <- survfit.o
}

par(mar=c(2,2,2,2))
# for clustering results using lym DMCTs:
plot(surv.res.l[[1]]$survfit, col = c("firebrick","dodgerblue","orange"),lwd=2, mark.time=TRUE, xlab="Years", ylab="OS",xscale = 365.25,cex.lab=0.6,cex.axis=0.6,mgp = c(1, 0.5, 0),tck = -0.03,main="Lym clusters",cex.main=0.7) 
text(x = 1400,y=0.1,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[1]," Chisq P=",signif(surv.res.l[[1]]$`pair-Pval`[1],2)),cex = 0.6)
text(x = 1400,y=0.15,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[2]," Chisq P=",signif(surv.res.l[[1]]$`pair-Pval`[2],2)),cex = 0.6)
text(x = 1400,y=0.2,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[3]," Chisq P=",signif(surv.res.l[[1]]$`pair-Pval`[3],2)),cex = 0.6)

# for clustering results using Hep DMCTs:
plot(surv.res.l[[2]]$survfit, col = c("#FF8080","#80FF80","#8080FF"),lwd=2, mark.time=TRUE, xlab="Years", ylab="OS",xscale = 365.25,cex.lab=0.6,cex.axis=0.6,mgp = c(1, 0.5, 0),tck = -0.03,main="Hep clusters",cex.main=0.7) 
text(x = 1400,y=0.1,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[1]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[1],2)),cex = 0.6)
text(x = 1400,y=0.15,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[2]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[2],2)),cex = 0.6)
text(x = 1400,y=0.2,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[3]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[3],2)),cex = 0.6)

# for clustering results using EC DMCTs:
plot(surv.res.l[[3]]$survfit, col =  c("#F06000","#408030"),lwd=2, mark.time=TRUE, xlab="Years", ylab="OS",xscale = 365.25,cex.lab=0.6,cex.axis=0.6,mgp = c(1, 0.5, 0),tck = -0.03,main="EC clusters",cex.main=0.7) 
text(x = 1400,y=0.1,label=paste0(names(surv.res.l[[3]]$`pair-Pval`)[1]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[1],2)),cex = 0.6)


## ----chunk8, eval=F,echo=T----------------------------------------------------
#  res.m <- GenResidualMat(LIHC_data$DNAm,estCTF.m = LIHC_estF,standardize = T,ncores = 40)
#  jive.results.l <- DoCELTYC(res.m,method = "jive",maxK = 3,dmct.lv = selDMCT,title = "cluster-test")
#  jive.clust.v <- jive.results.l$Lym[[3]]$consensusClass # extract the clustering results on JIVE-derived individual matrix for lymphocyte-specific DMCTs, when cluster number=3
#  

## ----sessionInfo, eval=T, echo=T----------------------------------------------
sessionInfo()

