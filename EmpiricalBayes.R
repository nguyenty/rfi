#install.packages("rjags")
library(rjags)
set.seed(20141020)
#############modelm - using point mass mixture prior for signals###############
modelm <- "
model{
# likelihood 
for (i in 1:length(y)){
y[i] ~ dnorm(beta[i], 1/sigma[i]^2)
beta[i] <- (1 - bin.beta[i])*norm.beta[i]
bin.beta[i] ~ dbern(pi.beta)
norm.beta[i] ~ dnorm(beta0, 1/sigma0^2)
sigma[i] ~ dunif(0, 100)
}

# prior distribution for the parameters ####
pi.beta ~ dbeta(5, 1)
beta0 ~ dnorm(0, 100)
sigma0 ~ dunif(0, 100)
}
"
# load("U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/ks/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_fit.RData")
# load("U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/ks/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")

load("Model1_fit.RData")
load("Model1_result.RData")

#hist(fit$coef[,6], nclass = 1000)
new.beta <- list()
for(i in 1:23)
{ 
  data <- list(y = fit$coef[,i])
  
  m0 <- proc.time()
  mm <- jags.model(textConnection(modelm), data,n.chains = 1) # mix point mass
  resm <- coda.samples(mm, c("beta","sigma","beta0","sigma0","pi.beta", 
                             "bin.beta")
                       , 2000) # mix point mass
  
  sigma.i <- apply(resm[[1]][,24563:36842], 2, mean)
  sigma.0 <- mean(resm[[1]][, "sigma0"])
  beta.0 <- mean(resm[[1]][,"beta0"])
  
  new.beta[[i]] <- fit$coef[,i]*(1/sigma.i^2)/(1/sigma.i^2 + 1/sigma.0^2) + 
    beta.0 * (1/sigma.0^2)/(1/sigma.i^2 + 1/sigma.0^2)
  
  
}

save(new.beta, file = "new.beta.RData")

source("quasiseq shrinkage functions.R")
myQL.fit

new.beta.matrix <- matrix()
for(i in 1:23){
  new.beta.matrix <- cbind(new.beta.matrix, new.beta[[i]])
}
dim(new.beta.matrix)
str(new.beta)
length(new.beta[[2]])

load("new.beta.RData")
str(new.beta)


str(fit)
#### 

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

#resultdir <- '/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired'
resultdir <- "U:/R/RA/Data/RFI-newdata/resultpairedlogcbc"
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
row.names(scount) <- scount[,1]
# dim(scount)
# str(scount)
# which(scount[,1] %in%"ENSSSCG00000007978")
# which(scount[,1] %in%"ENSSSCG00000014725")
# 
# scount[which(scount[,1] %in%"ENSSSCG00000007978"), ]
# scount[which(scount[,1] %in%"ENSSSCG00000014725"), ]

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
