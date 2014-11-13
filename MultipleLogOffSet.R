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


load("new.beta.RData")





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
list_model <- function(full_model){
  n <- dim(full_model)[2]
  variable_name <- colnames(full_model)[-1]
  variable_name <- gsub(":", "", variable_name)
  for (i in 2:8){variable_name <- gsub(i, "", variable_name)}
  test.mat <- NULL
  design.list <- vector("list", n)
  design.list[[1]] <- full_model
  for (i in 2:n) {
    design.list[[i]] <- as.matrix(full_model[,-i])
    test.mat <- rbind(test.mat, c(1,i))
  }
  
  row.names(test.mat) <-  variable_name
  
  if (any(variable_name == "Block") & any(variable_name == "Blockorder")){
    ind_block <- which(variable_name == "Block")[1]
    ind_blockorder <- which(variable_name == "Blockorder")[1] 
    nlist <- n - sum(variable_name == "Block") -
      sum(variable_name == "Blockorder") + 2
    design.list <- vector("list", nlist)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(nlist-1), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 1:(ind_block-1)){ # i <- ind_block - 1
      design.list[[i+1]] <- as.matrix(full_model[,-(i+1)])
      row.names(test.mat)[i] <- variable_name[i]
    }
    
    design.list[[ind_block + 1]] <- full_model[,-((ind_block+1):(ind_block+3))]
    
    design.list[[ind_block +2]] <- full_model[,-((ind_blockorder+1):(ind_blockorder+7))]
    row.names(test.mat)[ind_block + 1] <- "Blockorder" 
    if(ind_blockorder + 7 < n){ # i <- 15 # colnames(design.list[[i]])
      for(i in ((ind_block+3):(n-8))){
        design.list[[i]] <- full_model[,-(i + 8)] # colnames (full_model)
        row.names(test.mat)[i-1] <- variable_name[i+7]    
      }
    }
  }
  
  if (any(variable_name == "Block") & all(variable_name != "Blockorder")){
    ind_block <- which(variable_name == "Block")[1]
    nlist <- n - sum(variable_name == "Block") + 1
    design.list <- vector("list", nlist)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(nlist-1), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 1:(ind_block-1)){ # i <- ind_block - 1
      design.list[[i+1]] <- as.matrix(full_model[,-(i+1)])
      row.names(test.mat)[i] <- variable_name[i]
    }
    
    design.list[[ind_block + 1]] <- full_model[,-((ind_block+1):(ind_block+3))]
    if(ind_block + 3 < n){ # i <- 15 # colnames(design.list[[i]])
      for(i in ((ind_block+2):(n-2))){
        design.list[[i]] <- full_model[,-(i + 2)] # colnames (full_model)
        row.names(test.mat)[i-1] <- variable_name[i+1]    
      }
    }
  }
  
  if (all(variable_name != "Block") & any(variable_name == "Blockorder")){
    ind_blockorder <- which(variable_name == "Blockorder")[1] 
    nlist <- n - sum(variable_name == "Blockorder") + 1
    design.list <- vector("list", nlist)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(nlist-1), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 1:(ind_blockorder-1)){ # i <- ind_blockorder - 1
      design.list[[i+1]] <- as.matrix(full_model[,-(i+1)])
      row.names(test.mat)[i] <- variable_name[i]
    }
    
    design.list[[ind_blockorder + 1]] <- full_model[,-((ind_blockorder+1):(ind_blockorder+7))]
    row.names(test.mat)[ind_blockorder] <- "Blockorder" 
    if(ind_blockorder + 7 < n){ # i <- 15 # colnames(design.list[[i]])
      for(i in ((ind_blockorder+2):(n-6))){
        design.list[[i]] <- full_model[,-(i + 6)] # colnames (full_model)
        row.names(test.mat)[i-1] <- variable_name[i+5]    
      }
    }
  }
  if (n ==2) design.list[[2]] <- rep(1, nrow(full_model))
  return(list(design.list = design.list, test.mat = test.mat))
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

bvec <- fit$coef[-del.gene, 3]
svec <- ses[-del.gene,3]^2
length(bvec)
length(svec)

counts2 <- counts[-del.gene,]
log.offset <- log(apply(counts2, 2, quantile, 0.75))
model_th <- 111111 # test for Emperial Bayes method#####
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block + Blockorder)


list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
fit2 <- QL.fit(counts2, design.list, test.mat, # dim(counts)
               log.offset = log.offset, print.progress=TRUE,
               Model = "NegBin")

result2<- QL.results(fit2, Plot = FALSE)



ses <- matrix(0, nrow = nrow(fit2$coefficients), ncol = ncol(fit2$coefficients))

for(j in 1:nrow(fit2$coefficients)){ 
  se.new <- se.b1(fit2$coefficients[j,],fit2$NB.disp[j],full_model)
  ses[j,] <- se.new
}
dim(ses)

del.gene <- which(apply(ses, 1, sum)==0)

bvec <- fit$coef[-del.gene, 3]
svec <- ses[-del.gene,3]^2
length(bvec)
length(svec)
#### fit QL.fit model excluding those genes with non inverse Hessian Matrix#####

## first: obtain estimate of sigma_k for each sample, k = 1, 23#####

s0b0 <- function(x){
  s0 <- x[1]; b0 <- x[2]
  l <- sum(-log(svec+s0) - (bvec - b0)^2/(svec+s0))/2
  -l
}
g.s0b0 <- function(x){
  s0 <- x[1]; b0 <- x[2]
  c(sum(-.5/(svec + s0)+.5/(svec + s0)^2*(bvec-b0)^2), 
  sum(-1/(svec + s0)^2*(bvec-b0)))
}
# k <- 1 # s0 <- 1; b0 <- 0

optimx(c(1, 0), s0b0, g.s0b0, method = "BFGS")
# 
# optimx(c(1, 0), s0b0)
# install.packages("optimx")
# library("optimx")




load("new.beta.RData")
str(new.beta)

bvec <- matrix(0, nrow = length(new.beta[[1]]), ncol = length(new.beta))

for(i in 1:length(new.beta)){
  bvec[,i] <- new.beta[[i]]
}

model_th <- 11111 # test for Emperial Bayes method#####
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block + Blockorder)
new.offset <- t(full_model[, -2]%*% t(bvec[,-2]))
dim(new.offset)

source("quasiseq shrinkage functions2.R")


design.list <- vector("list", 2)
design.list[[1]] <- model.matrix(~Line)
design.list[[2]] <- rep(1, 31)


log.offset <- log(apply(counts, 2, quantile, .75))
fit <- myQL.fit(counts, design.list, # dim(counts)
              log.offset = log.offset, betavec = new.offset, print.progress=TRUE,
              Model = "NegBin", method = "optim")

result<- QL.results(fit, Plot = FALSE)
str(result)
hist(result$P.values[[3]])
sum(result$Q.values[[3]]<=0.05)
