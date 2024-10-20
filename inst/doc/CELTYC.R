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
#  library(EpiDISH)
#  phe.v <- pheno.df$cancer # a vector indicating whether a sample is "Cancer" or "Normal"
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
#  celldmc.o <- CellDMC(beta.m = bmiq.m[,-na.idx], pheno.v=phe.v[-na.idx], frac.m  = estF.m[-na.idx,], mc.cores = 40,cov.mod = cov.mod.use)
#  dmct.lv <- GetDMCT(celldmc.o,nDMCT = 100,specific = T,nSpDMCT = 100,common = T,nCT = 3,nCmDMCT = 100) # nSpDMCT is the threshold to save DMCTs specific to a particular cell type, and nCmDMCT is the threshold to save DMCTs common to nCT cell types

## ----chunk4, message=F--------------------------------------------------------
library(SuperExactTest)
supertest.o <- supertest(LIHC_data$allDMCT$dmct,n=403197)
plot(supertest.o,"landscape",sort.by = "size",color.scale.cex = 0.8,overlap.size.cex = 0.7,cex=0.8)

## ----chunk5,echo=T------------------------------------------------------------
res.m <- GenResidualMat(LIHC_data$DNAm,estCTF.m = LIHC_data$estF,standardize = T,ncores = 40) # regress out CTFs from DNAm matrix and get standardized residual matrix
CELTYC.results.l <- DoCELTYC(res.m,method = "consensus",maxK = 3,dmct.lv = LIHC_data$allDMCT$spec[c("Lym","Hep","EC")],title = "cluster-test") # do consensus clustering on standardized residual matrix over DMCTs specific to lymphocytes, hepatocytes and endothelial cells

## ----chunk6,message=F---------------------------------------------------------
library(ComplexHeatmap)
lym.clust.v <- CELTYC.results.l$Lym[[3]]$consensusClass
hep.clust.v <- CELTYC.results.l$Hep[[3]]$consensusClass
ec.clust.v <- CELTYC.results.l$EC[[2]]$consensusClass # if the chosen cluster number is 3, the 3rd cluster contains few samples
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
  # plot Kaplan Meier curves 
  surv.o <- Surv(time = LIHC_data$pheno$OS.time,event = LIHC_data$pheno$event)
  survfit.o <- survfit(surv.o ~ clust.v) 
  surv.res.l[[m]] <- list()
  surv.res.l[[m]][["pair-Pval"]] <- chi.surv.p.v
  surv.res.l[[m]][["survfit"]] <- survfit.o
}

par(mar=c(2,2,2,2))
# for clustering results using lym DMCTs:
plot(surv.res.l[[1]]$survfit, col = c("firebrick","dodgerblue","orange"),lwd=2, mark.time=TRUE, xlab="Years", ylab="OS",xscale = 365.25,cex.lab=0.8,cex.axis=0.8,mgp = c(1, 0.5, 0),tck = -0.03,main="Lym clusters",cex.main=1) 
text(x = 1400,y=0.1,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[1]," Chisq P=",signif(surv.res.l[[1]]$`pair-Pval`[1],2)),cex = 0.9)
text(x = 1400,y=0.15,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[2]," Chisq P=",signif(surv.res.l[[1]]$`pair-Pval`[2],2)),cex = 0.9)
text(x = 1400,y=0.2,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[3]," Chisq P=",signif(surv.res.l[[1]]$`pair-Pval`[3],2)),cex = 0.9)
legend("topright",legend=paste0("cl",1:3),col=c("firebrick","dodgerblue","orange"),
       lty = 1,inset = 0.05,cex = 0.9,bty="n")


# for clustering results using Hep DMCTs:
plot(surv.res.l[[2]]$survfit, col = c("#FF8080","#80FF80","#8080FF"),lwd=2, mark.time=TRUE, xlab="Years", ylab="OS",xscale = 365.25,cex.lab=0.8,cex.axis=0.8,mgp = c(1, 0.5, 0),tck = -0.03,main="Hep clusters",cex.main=1) 
text(x = 1400,y=0.1,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[1]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[1],2)),cex = 0.9)
text(x = 1400,y=0.15,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[2]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[2],2)),cex = 0.9)
text(x = 1400,y=0.2,label=paste0(names(surv.res.l[[1]]$`pair-Pval`)[3]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[3],2)),cex = 0.9)
legend("topright",legend=paste0("cl",1:3),col=c("#FF8080","#80FF80","#8080FF"),
       lty = 1,inset = 0.05,cex = 0.9,bty="n")

# for clustering results using EC DMCTs:
plot(surv.res.l[[3]]$survfit, col =  c("#F06000","#408030"),lwd=2, mark.time=TRUE, xlab="Years", ylab="OS",xscale = 365.25,cex.lab=0.8,cex.axis=0.8,mgp = c(1, 0.5, 0),tck = -0.03,main="EC clusters",cex.main=1) 
text(x = 1400,y=0.1,label=paste0(names(surv.res.l[[3]]$`pair-Pval`)[1]," Chisq P=",signif(surv.res.l[[2]]$`pair-Pval`[1],2)),cex = 0.9)
legend("topright",legend=paste0("cl",1:2),col=c("#F06000","#408030"),
       lty = 1,inset = 0.05,cex = 0.9,bty="n")


## ----chunk8, eval=F,echo=T----------------------------------------------------
#  selDMCT.lv <- c(LIHC_data$allDMCT$spec[c("Lym","Hep","EC")],LIHC_data$allDMCT$comm["EC_Hep_Lym"]) # extract the cell-type specific DMCTs and the common DMCTs
#  jive.results.l <- DoCELTYC(data.m = res.m,method = "jive",maxK = 3,dmct.lv = selDMCT.lv,title = "cluster-test")
#  jive.IV.lym.clust.l <- jive.results.l$Lym # extract the clustering results on JIVE-derived individual variation matrix for lymphocyte-specific DMCTs

## ----chunk9, eval=T,echo=T----------------------------------------------------
jive.clust.v <- LIHC_data$jive_IV_lym[[2]]$consensusClass # extract the cluster assignments when cluster number is 2
surv.o <- Surv(time = LIHC_data$pheno$OS.time,event = LIHC_data$pheno$event)
cox.o <- coxph(surv.o~jive.clust.v)
chisq.pval <- pchisq(cox.o$score,1,lower.tail = F)
survfit.o <- survfit(surv.o ~ jive.clust.v)
par(mar=c(2,2,2,2))
plot(survfit.o, col =  c("coral","#CF30CF"),lwd=2, mark.time=TRUE, xlab="Years", ylab="OS",xscale = 365.25,cex.lab=0.8,cex.axis=0.8,mgp = c(1, 0.5, 0),tck = -0.03,main="JIVE IV Lym clusters",cex.main=1)
text(x = 1400,y=0.1,label=paste0("cl 1,2 Chisq P=",signif(chisq.pval,2)),cex = 0.9)
legend("topright",legend=paste0("cl",1:2),col=c("coral","#CF30CF"),
       lty = 1,inset = 0.05,cex = 0.9,bty="n")

## ----sessionInfo, eval=T, echo=T----------------------------------------------
sessionInfo()

