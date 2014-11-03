require(Matrix)
#library(QuasiSeq)
library(edgeR)
require(reshape)
require(plyr)
library(fields)
library(reshape)
library(fdrtool)
source("QL.fit.R")
source("NBDev.R")
source("PoisDev.R")
source("QL.results.R")
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
row.names(scount) <- scount[,1]
scount <- scount[-c(which(scount[,1] %in%"ENSSSCG00000007978"),
                    which(scount[,1] %in%"ENSSSCG00000014725")),]


cbc <- read.table('CBC data for pigs with RNA-seq data avaible.txt',
                  header =T)

metadata <- read.table("meta_data_RNA-seq_G9P2.txt", 
                       header = T)

rfiadj <- read.csv("g8p2_g9p2-rfiadj-FINAL_Jan_2014_rfiadjusted.csv", 
                   header = T)

##### cleaning data####



cbc <- cbc[order(cbc$ear), ]
metadata <- metadata[order(metadata$idpig), ]
rfiadj <- rfiadj[order(rfiadj$idpig),]

fullidpig <- as.numeric(paste("20900", metadata$idpig, sep = ""))

covset <- cbind(metadata[, -4], rfiadj[rfiadj$idpig %in% fullidpig, c("rfi.ADJUSTED")], 
                cbc[, c("iddam", "idsire", "Neutrophil", "Lymphocyte", "Monocyte",
                        "Eosinophil", "Basophil" )])
colnames(covset) <- c("idpig", "Line", "Diet",  "Block", "Blockorder", "Concb", 
                      "RINb", "Conca", "RINa", "RFI",
                      "iddam", "idsire", "neut",   
                      "lymp", "mono","eosi", "baso")

#####set of covariates considered ####
covset <- cbind(covset, lneut = log(covset$neut), llymp = log(covset$lymp), 
                lmono = log(covset$mono), leosi = log(covset$eosi), 
                lbaso = log(covset$baso)) 
covset$Line <- as.factor(covset$Line)
covset$Diet <- as.factor(covset$Diet)
covset$Block <- as.factor(covset$Block)
covset$Blockorder <- as.factor(covset$Blockorder)
covset$iddam <- as.factor(covset$iddam)
covset$idsire <- as.factor(covset$idsire)
levels(covset$idsire) <- 1:11
levels(covset$iddam) <- 1:20

covset[, c("iddam", "idsire")]
#detach(covset)
attach(covset)

# counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
#                              rowMeans(scount[,-1])>8 & 
#                              rowSums(scount[,-1][,Line ==1] > 0) >0 &
#                              rowSums(scount[,-1][, Line ==2] >0) >0 ,-1])

counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>8 ,-1])

load("U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/pvalue05/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_fit.RData")
load("U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/pvalue05/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")

dim(fit$coef)
#dl2: design MATRIX; b: fit coefficients, w: fit NB.disp
f <- function(m) class(try(solve(m),silent=T))=="matrix"
se.b1 <- function(b,w,dl2){
  
  X <- dl2
  
  eta <- c(X%*%b)
  
  mu <- exp(eta)
  
  fish <- t(X)%*%(diag(nrow(X))*mu/(mu*w+1))%*%X
  
  if(f(fish)){
    fish.inv <- solve(fish)
    return(sqrt(diag(fish.inv)))
  }else{
    return(rep(0, ncol(X)))
  }
}


prefit <- fit
model_th <- 1
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block + Blockorder)

#prefit is the output of QL.fit using the full model including covariates
# 
# b <- prefit$coefficients[j,]
# w <- prefit$NB.disp[j]
# dl2 <- full_model

ses <- matrix(0, nrow = nrow(prefit$coefficients), ncol = ncol(prefit$coefficients))

for(j in 1:nrow(prefit$coefficients)){ 
  se.new <- se.b1(prefit$coefficients[j,],prefit$NB.disp[j],full_model)
  ses[j,] <- se.new
}
dim(ses)

del.gene <- which(apply(ses, 1, sum)==0)


j <- 2647 # 7299, 9257, 11574
j <- 1
se.b1(prefit$coefficients[j,],prefit$NB.disp[j],full_model)
 load("new.beta.RData")
str(new.beta[[1]])
dim(full_model)
colnames(full_model)
full_model[, -c(1, 2)]
#### fit QL.fit model excluding those genes with non inverse Hessian Matrix#####

## first: obtain estimate of sigma_k for each sample, k = 1, 23#####

s0b0 <- function(x){
  s0 <- x[1]; b0 <- x[2]
  l <- sum(-log(svec+s0) - (bvec - b0)^2/(svec+s0))
  -l
}
g.s0b0 <- function(x){
  s0 <- x[1]; b0 <- x[2]
  c(sum(1/(svec + s0)+1/(svec + s0)^2*(bvec-b0)^2), 
  sum(2/(svec + s0)^2*(bvec-b0)^2))
}
# k <- 1 # s0 <- 1; b0 <- 0
bvec <- fit$coef[-del.gene, 1]
svec <- ses[-del.gene,1]
length(bvec)
length(svec)
optim(c(1, 0), s0b0, g.s0b0, method = "BFGS")
?optim
