library(Matrix)
# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(edgeR)
library(reshape)
library(plyr)
library(fields)
library(reshape)
library(fdrtool)
library(QuasiSeq)

#resultdir <- '/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired'
#resultdir <- "/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsimulation"
resultdir <- "U:/R/RA/Data/RFI-newdata/resultsimulation"
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
scount <- scount[-c(which(scount[,1] %in%"ENSSSCG00000007978"),
                    which(scount[,1] %in%"ENSSSCG00000014725")),]
counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>8 ,-1])


###List of models function ####
covset <- read.table("covset.txt")
attach(covset)
Blockorder <- as.factor(Blockorder)
Block <- as.factor(Block)
Line <- as.factor(Line)
Diet <- as.factor(Diet)
load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue052/Model777.Line.Concb.RINa.neut.lymp.mono.baso.Block/Model777_fit.RData")
load("U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue052/Model777.Line.Concb.RINa.neut.lymp.mono.baso.Block/Model777_result.RData")
# 
# load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpairedlogcbc/pvalue05/Model7.Line.Concb.RINa.lneut.llymp.lmono.lbaso.Block/Model7_fit.RData")
# load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpairedlogcbc/pvalue05/Model7.Line.Concb.RINa.lneut.llymp.lmono.lbaso.Block/Model7_result.RData")

full_model <- model.matrix(~Line + Concb + RINa + 
                             lneut + llymp + lmono + 
                             lbaso + Block)
coef_beta <- fit$coef 
ee_coef <- sort(result$Q.values[[3]][,"Line"], index.return = T)$ix[(dim(coef_beta)[1]-result$m0[3,1]+1):dim(coef_beta)[1]]
coef_beta[ee_coef,2] <- 0
cor_fit_count <- laply(1:dim(coef_beta)[1], function(i)
  cor(fit$fitted.values[i,], counts[i,]))
used_gene <- which(cor_fit_count >.8)#length(used_gene)
used_beta <- coef_beta[used_gene,]
used_omega <- fit$NB.disp[used_gene]
used_count <- counts[used_gene,]
degene <- which(used_beta[, 2] !=0) # length(degene)
mu <- fit$fitted[used_gene,]
# degene
# length(used_gene)
