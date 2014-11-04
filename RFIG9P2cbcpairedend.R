require(Matrix)
library(QuasiSeq)
library(edgeR)
require(reshape)
require(plyr)
library(fields)
library(reshape)
library(fdrtool)
# source("QL.fit.R")
# source("NBDev.R")
# source("PoisDev.R")
# source("QL.results.R")

#resultdir <- '/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired'
resultdir <- "U:/R/RA/Data/RFI-newdata/resultpairedcbc"
scount <- read.table("paired end uniquely mapped reads count table.txt", 
                     header = T)
row.names(scount) <- scount[,1]
# dim(scount)
# scount
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

write.csv(counts, file = "counts_12280.csv")
dim(counts)
dim(scount)
dim(counts)
log.offset <- log(apply(counts, 2, quantile, .75))
g_cdf <- function(z){
  e <- ecdf(z)
  g <- grenander(e)
  g
}

sel_criteria <- function(result){
  dat <- result$P.values[[3]][,colnames(result$P.values[[3]])]
  # Crames Von Miser statistics
  dat <- as.matrix(dat)
  cvm <- apply(dat, 2, function(z)sum((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2 *
                                        diff(c(0,g_cdf(z)$x.knots))))
  # Kolmogorow Smirnov statistics 
  ks <- apply(dat, 2, function(z)max((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2))
  
  # Anderson-Darling statistics
  ad <- apply(dat, 2, function(z)sum((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2/
                                       g_cdf(z)$x.knots*(1-g_cdf(z)$x.knots)*
                                       diff(c(0,g_cdf(z)$x.knots))))
  # Proportion of pvalue less than 0.05
  pvalue_05 <- apply(dat<=0.05, 2, sum)
  
  out <- data.frame(pvalue05 = order(pvalue_05),
                    ad = order(ad),
                    cvm = order(cvm),
                    ks = order(ks))
  
  return(out)
}


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

#### clean list_model function  #####

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


fit_model <- function(full_model, model_th, criteria){ # model_th <- 1
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, # dim(counts)
                log.offset = log.offset, print.progress=TRUE,
                Model = "NegBin", method = "optim")
  result<- QL.results(fit, Plot = FALSE)
  res_sel <- sel_criteria(result)
  k <- nrow(test.mat)
  name_model <- NULL 
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(resultdir,"/", colnames(res_sel)[criteria], "/Model",model_th,name_model, sep ="")
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
  
  return(res_sel)
}


out_pairedend_cbc <-  data.frame(Date=as.Date(character()),
                                   File=character(), 
                                   User=character(), 
                                   stringsAsFactors=FALSE) 

for(i in 1){ # i <- 1
  model_th <- 1
  full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                               neut + lymp + mono + eosi + baso + 
                               Block + Blockorder)
  repeat{
    pm1 <- proc.time()
    out_model <- fit_model(full_model, model_th, i)
    assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
    proc.time() -pm1
    ms_val <- get(paste("ms_criteria", model_th, sep = "_" ))
    cov_del <- ms_val[1,i] # cov_del <- 14; i <- 1
    
    cov_set <- list_model(full_model)$test.mat # dim(cov_set)
    res <- data.frame(criteria = colnames(ms_val)[i], model = model_th, cov_del = rownames(cov_set)[ cov_del])
    out_pairedend_cbc <- rbind(out_pairedend_cbc, res)
    if (cov_del ==1) break
    block_ind <- grep("Block2", colnames(full_model))
    blockorder_ind <-grep("Blockorder", colnames(full_model))
    
    if(length(block_ind)!=0){
      if(cov_del+1 == cov_set["Block",2])
        full_model <- full_model[, -c(block_ind, block_ind+1, block_ind+2)]
    }
    
    if(length(blockorder_ind)!=0){
      if(cov_del+1 == cov_set["Blockorder", 2])
        full_model <- full_model[, -blockorder_ind]
    }
    
    full_model <- full_model[, -(cov_del +1)]
    model_th <- model_th +1
  }
}

write.csv(out_pairedend_cbc, "out_pairedend_cbc.csv", row.names = F)



model_th <- 7777 ## Using QuasiSeq1.0.4#####
full_model <- model.matrix(~Line + Concb   + RINa + 
                             neut + lymp + mono  + baso + 
                             Block )

out_model <- fit_model(full_model, model_th, 1)
## rerun above model 777 with new QuasiSeq####
## has some error on glm.fit ###

list_out <- list_model(full_model)[[1]][c(1, 2)]

design.list <- list_out

fit2 <- QL.fit(counts, design.list,  # dim(counts)
              log.offset = log.offset, print.progress=TRUE,
              Model = "NegBin", method = "optim")
result2<- QL.results(fit2, Plot = FALSE)
(result[[3]])
str(result)
sum((result2$Q.values[[3]]<=.05)&(abs(fit2$coef[,2]/log(2))>=1))
sum(result2$Q.values[[3]]<=.05)
hist(result2$P.values[[3]], nclass  =100)
