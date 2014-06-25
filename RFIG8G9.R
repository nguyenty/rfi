require(Matrix)
#library(QuasiSeq)
library(edgeR)
require(reshape)
require(plyr)
library(fields)
library(reshape)
dir.source <- "U:/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
#dir.source <- "/home/ntyet/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
source(paste(dir.source, "QL.fit.R",sep=""))
source(paste(dir.source, "NBDev.R",sep =""))
source(paste(dir.source, "PoisDev.R",sep =""))
source(paste(dir.source, "QL.results.R",sep =""))

### Reading data #######
# scount <- read.table("/home/ntyet/research/RFI-newdata/Data for Yet/single end uniquely mapped reads count table for Yet.txt", 
#                      header = T)
# cbc <- read.table('/home/ntyet/research/RFI-newdata/Data for Yet/CBC data for pigs with RNA-seq data avaible.txt',
#                   header =T)
# 
# metadata <- read.table("/home/ntyet/research/RFI-newdata/Data for Yet/meta_data_RNA-seq_G9P2.txt", 
#                        header = T)
# 
# rfiadj <- read.csv("/home/ntyet/research/RFI-newdata/Data for Yet/g8p2_g9p2-rfiadj-FINAL_Jan_2014_rfiadjusted.csv", 
#                      header = T)

#resultdir <- "/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result"
resultdir <- "U:/R/RA/Data/RFI-newdata/result"
scount <- read.table("single end uniquely mapped reads count table for Yet.txt", 
                     header = T)
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

###List of models function ####
### Case 1: no cbc data ####
# dim(covset)
# colnames(covset)
# full_model <- model.matrix(~Line*Diet*RFI + Block+ Blockorder)
# colnames(full_model)
# rankMatrix(full_model)
# list_model(full_model)
log.gamma <- function(counts, disp){
  log.g <- NULL
  n <- length(counts)
  for(i in 1:n){
    log.g[i] <- sum(log(0:(counts[i]-1)+disp)) - sum(log(1:counts[i]) )
  }
  return(log.g)  
}

SAT.LIKE2<-function(count,disp){
  means<-count
  like<-disp*log(disp/(disp+means))
  like[count!=0]<-like[count!=0]+count[count!=0]*log(means[count!=0]/(disp+means[count!=0]))+
    log.gamma(count[count!=0],disp )
  sum(like)
}

## Function to calculate AIC of the QL.fit model 

AIC.QL <- function(counts,QL.fit.object){
  n <- dim(counts)[2]
  m <- dim(counts)[1]
  disp <- 1/QL.fit.object$NB.disp
  den.df <- QL.fit.object$den.df
  phi.hat.dev <- QL.fit.object$phi.hat.dev
  p <- n - den.df
  dev <- phi.hat.dev*den.df
  L0 <- NULL
  for (i in 1:m){
    L0[i] <- SAT.LIKE2(counts[i,],disp[i])
  }
  
  return(dev-2*L0+2*p)
}

counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                           rowMeans(scount[,-1])>1 & 
                             rowSums(scount[,-1][,Line ==1] > 0) >0 &
                             rowSums(scount[,-1][, Line ==2] >0) >0 ,-1])


dim(scount)
dim(counts)
log.offset <- log(apply(counts, 2, quantile, .75))
###############

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
    # row.names(test.mat)
    #colnames(design.list[[ind_block+1]])
    design.list[[ind_block +2]] <- full_model[,-((ind_blockorder+1):(ind_blockorder+7))]
    row.names(test.mat)[ind_block + 1] <- "Blockorder" 
    # colnames(design.list[[ind_block+2]])
    # colnames(design.list[[1]])
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
    # row.names(test.mat)[ind_block] <- "Block"
    # row.names(test.mat)
    #colnames(design.list[[ind_block+1]])
#     design.list[[ind_block +2]] <- full_model[,-((ind_blockorder+1):(ind_blockorder+7))]
#     row.names(test.mat)[ind_block + 1] <- "Blockorder" 
    # colnames(design.list[[ind_block+2]])
    # colnames(design.list[[1]])
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
  # row.names(test.mat)[[ind_blockorder]]
  # colnames(design.list[[ind_blockorder+1]])
  row.names(test.mat)[ind_blockorder] <- "Blockorder" 
  # colnames(design.list[[ind_block+2]])
  # colnames(design.list[[1]])
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


## Function do all the things with input Full model

fit_model <- function(full_model, model_th){ # model_th <- 1
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, 
                log.offset = log.offset, print.progress=FALSE,
                Model = "NegBin")
  result<- QL.results(fit, Plot = FALSE)
  k <- nrow(test.mat)
  name_model <- NULL
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(resultdir, "/Model",model_th,name_model, sep ="")
  dir.create(model_dir, showWarnings = FALSE)
  save(result, file = paste(model_dir,"/Model",model_th, "_result.RData", sep =""))
  save(fit, file = paste(model_dir,"/Model",model_th, "_fit.RData", sep =""))
  for(i in 1:(nrow(test.mat))){
    postscript(paste(model_dir,"/Model", 
                     model_th, row.names(test.mat)[i],".eps", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
    
    pdf(paste(model_dir,"/Model", 
              model_th, row.names(test.mat)[i],".pdf", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
  }
  print(paste("Model", model_th, sep = " "))
  
  return(list(mean_model = mean(fit$phi.hat.dev), 
              median_model = median(fit$phi.hat.dev), 
              AIC_model = mean(AIC.QL(counts, fit))))
}

# Model 1
m <- 1
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block + Blockorder)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 2
m <- 2
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + Concb + 
                             RINb + Conca + RINa + 
                             Block + Blockorder)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 




# Model 3
m <- 3
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + Concb + 
                             RINb + Conca + RINa + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 4
m <- 4
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + Concb + 
                             RINb + RINa + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 5
m <- 5
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + Concb + 
                              RINa + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 




# Model 6
m <- 6
model_th <- m
full_model <- model.matrix(~Line*Diet + Diet*RFI + Line*RFI + Concb + 
                             RINa + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 7
m <- 7
model_th <- m
full_model <- model.matrix(~Line + Diet + Diet*RFI + Line*RFI + Concb + 
                             RINa + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 8
m <- 8
model_th <- m
full_model <- model.matrix(~Line + Diet + Diet*RFI + Line*RFI + Concb + 
                             RINa + 
                             Block)[, -3]
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 




# Model 9
m <- 9
model_th <- m
full_model <- model.matrix(~Line + Diet + Diet*RFI + Line*RFI + Concb + 
                             RINa + 
                             Block)[, -c(3,4)]
#colnames(full_model)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 10
m <- 10
model_th <- m
full_model <- model.matrix(~Line + Diet + Diet*RFI + Line*RFI + Concb + 
                             RINa + 
                             Block)[, -c(3,4,11)]
#colnames(full_model)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 11
m <- 11
model_th <- m
full_model <- model.matrix(~Line + Concb + 
                             RINa + 
                             Block)
#colnames(full_model)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 100
m <- 100
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + 
                             Block + Blockorder)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 101
m <- 101
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + 
                             lneut + 
                             Block + Blockorder)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 

