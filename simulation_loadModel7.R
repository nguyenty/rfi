# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
library(Matrix)
library(edgeR)
library(reshape)
library(plyr)
library(fields)
library(reshape)
library(fdrtool)
library(QuasiSeq)
library(AUC)
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
scount <- scount[-c(which(scount[,1] %in%"ENSSSCG00000007978"),
                    which(scount[,1] %in%"ENSSSCG00000014725")),]
counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>8 ,-1])


###List of models function ####
covset <- read.table("covset.txt")
Blockorder <- as.factor(covset$Blockorder)
Block <- as.factor(covset$Block)
Line <- as.factor(covset$Line)
Diet <- as.factor(covset$Diet)
RFI <- covset$RFI
RINa <- covset$RINa
RINb <- covset$RINb
Conca <- covset$Conca
Concb <- covset$Concb
neut <- covset$neut
lymp <- covset$lymp
mono <- covset$mono
baso <- covset$baso
eosi <- covset$eosi
load("Model7_fit.RData")
load("Model7_result.RData")

full_model <- model.matrix(~Line + Concb + RINa + 
                             neut + lymp + mono + 
                             baso + Block)

coef_beta <- fit$coef 
ee_beta <- sort(result$Q.values[[3]][,"Line"], index.return = T)$ix[(dim(coef_beta)[1]-result$m0["QLSpline","Line"]+1):dim(coef_beta)[1]]
coef_beta[ee_beta,2] <- 0
Xb <- coef_beta%*%t(full_model) # dim(Xb)
log.offset <- t(as.matrix(log(apply(counts, 2, quantile, 0.75)))) # dim(ox)
ox <- log.offset[rep(1:nrow(log.offset), times = dim(fit$fitted)[1]), ]
mu <- exp(Xb+ox)
omega <- fit$NB.disp
degene <- which(coef_beta[, 2] !=0) # length(degene)
#length(degene)/length(omega)