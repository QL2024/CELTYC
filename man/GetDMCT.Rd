% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetDMCT.R
\name{GetDMCT}
\alias{GetDMCT}
\title{A function that returns a list of disease-associated DMCs altered in specific cell types (i.e. DMCTs)}
\usage{
GetDMCT(
  beta.m,
  pheno.v,
  ctf.m,
  adjPMethod = "fdr",
  adjPThresh = 0.05,
  cov.mod = NULL,
  ncores = 4
)
}
\arguments{
\item{beta.m}{A beta value matrix with rows labeling the CpGs and columns labeling samples.}

\item{pheno.v}{A vector of phenotype. Both binary and continuous/oderinal phenotypes are allowed, but NA is not allowed in pheno.v.}

\item{ctf.m}{A matrix contains fractions of each cell-type. Each row labels a sample, with the same order of the columns in beta.m.
Each column labels a cell-type. Column names, which are the names of cell-types, are required. The rowSums of ctf.m should be 1 or close to 1.}

\item{adjPMethod}{The method used to adjust p values. The method can be any of method accepted by \code{\link{[stats]p.adjust}}.}

\item{adjPThresh}{A numeric value, default as 0.05. This is used to call DMCTs. For each cell-type respectively, the CpG with the adjusted p values less than this threshold will be reported as DMCTs (-1 or 1) in the 'dmct' matrix in the returned list.}

\item{cov.mod}{A design matrix from `model.matrix`, which contains other covariates to be adjusted. For example, input model.matrix(~ geneder, data = pheno.df) to adjust gender. Do not put cell-type fraction here!}

\item{ncores}{The number of cores to use, i.e. at most how many threads will run simultaneously.}
}
\value{
A list containing 2 lists. The 1st list is the returned results of function \code{\link[EpiDISH]{CellDMC}}.
The 2nd list is a list of identified DMCTs for each cell type.
}
\description{
This function employs a previously published algorithm CellDMC (https://doi.org/10.1038/s41592-018-0213-x)
}
