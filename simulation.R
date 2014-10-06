require(Matrix)

# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
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
# check if the condition of count avarage is sastified
# change to simulation setup
# estimate the number 
# the following function is the same as function pval.hist
# in the paper of Pounds et al. 2012, EBT, with the estimation of 
# density obtained from the paper by Korbinian Strimmer 2008
# 

pval.hist.grenander <- function(p.value){
  grenander.out <- grenander(ecdf(p.value))
  p.brks <- c(0, grenander.out$x.knots)
  b.edf <- c(0, grenander.out$F.knots)
  p.diffs <- diff(p.brks)
  h.cdf <- approx(p.brks, b.edf, xout = p.value)$y  # get the histogram-based CDF estimates from each p-value
  p.hist <- exp(log(diff(b.edf))-log(diff(p.brks))) # get the hight for each histogram bar
  pi0.hat <- min(p.hist)                            # get the pi0 estimate from histogram bar
  h.ebp <- approx(p.brks, pi0.hat/c(p.hist, p.hist[length(p.hist)]), xout = p.value)$y # get the conservative EBP interpolation 
  h.fdr <- exp(log(pi0.hat) + log(p.value) - log(h.cdf))                                     # Get the histogram based FDR estimate
  h.ebp[p.value==0] <- 0
  h.fdr[p.value==0] <- 0
  return(list( p.value = p.value,          # input p-value,
               h.cdf = h.cdf,              # the histogram Grenander based cdf estimate
               h.fdr = h.fdr,              # the histogram Grenander based FDR
               h.ebp = h.ebp,              # the histogram Grenander based EBP
               p.brks = p.brks,            # the p-value break-points of the histogram
               p.hist = p.hist,            # the heights of each histogram bar
               edf.brks = b.edf,           # the breaks points in the EDF of the histogram estimator
               pi0.hat = pi0.hat))         # the histogram Grenander based estimate of the proportion of tests with a true null hypothesis
}
load("Model7_fitdat2.RData")
load("Model7_resultdat2.RData")

pvalue_line <- result$P.values[[3]][,"Line"]
qvalue_line <- result$Q.values[[3]][,"Line"]
gre_out <- pval.hist.grenander(pvalue_line)
ebp_line <- gre_out$h.ebp
fdr_line <- gre_out$h.fdr

full_model <- model.matrix(~Line + Concb + RINa + lneut + llymp + lmono + lbaso + Block)
coef_beta <- fit$coef 
coef_beta[,2] <- fit$coef[,2]*(ebp_line<0.5)

cor_fit_count <- laply(1:dim(coef_beta)[1], function(i)
  cor(fit$fitted.values[i,], counts[i,]))
set.seed(1)
used_gene <- which(cor_fit_count >.8)
used_beta <- coef_beta[used_gene,]
J <- dim(used_beta)[1]
s <- sample(J, J)
s <- s[order(s)]
used_omega <- fit$NB.disp[used_gene]
used_count <- counts[used_gene,]
degene <- which(ebp_line[used_gene]<0.5)
mu <- fit$fitted[used_gene,]
s_mu <- mu[s,]
s_omega <- used_omega[s]
s_degene <- intersect(degene, s) 

##Sim counts data
y <- array(0, dim = c(J,31))

save(sim_output, file = "sim_output.RData")

for(j in 1:J){
  repeat{
      for(k in 1:31){
        y[j,k]  <- rnbinom(n=1, size=1/s_omega[j], mu=s_mu[j,k])
      }
      if (mean(y[j,])>8& sum(y[j,]>0)>3) break
    }
}

sim_output <- list(used_gene = used_gene, used_omega = used_omega, used_count = used_count, 
                   degene_original = degene, s_degene = s_degene, y = y)



g_cdf <- function(z){
  e <- ecdf(z)
  g <- grenander(e)
  g
}
sel_criteria <- function(result){
  dat <- result$P.values[[3]][,colnames(result$P.values[[3]])]
  # Crames Von Miser statistics
  if(is.vector(dat)) dat <- as.matrix(dat)
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


## Function do all the things with input Full model

fit_model <- function(full_model, model_th, criteria){ # model_th <- 1
  counts <- y
  log.offset <- log(apply(y, 2, quantile, 0.75))
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, # dim(counts)
                log.offset = log.offset, print.progress=TRUE,
                Model = "NegBin")
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




list_cov_out <- df <- data.frame(Date=as.Date(character()),
                                 File=character(), 
                                 User=character(), 
                                 stringsAsFactors=FALSE) 

for(i in 1:4){ # i <- 2
  model_th <- m
  full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                               lneut + llymp + lmono + leosi + lbaso + 
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
    list_cov_out <- rbind(list_cov_out, res)
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
list_cov_out

                     
load("Model7_fitdat2.RData")
load("Model7_resultdat2.RData")
                     
str(fit)
dev_real <- fit$phi.hat.pearson
hist(dev_real[used_gene], nclass = 100)
summary(dev_real)
NB_real <- fit$NB.disp
hist(NB_real, nclass = 100)
summary(dev_real)
plot(y= dev_real, x= NB_real)
LRT<- fit$LRT[, 1]
plot(1-pchisq(LRT, df=1), result$P.values[[3]][, 1], xlim = c(0, 0.05), 
     ylim = c(0, 0.05))
lines(c(0, 1), c(0, 1), col = 2)
sum(1-pchisq(LRT, df=1)<=.01)
sum(result$P.values[[3]][, 1] <= .01)
hist(1-pchisq(LRT, df=1))
hist(result$P.values[[3]][, 1])
pvalue <- function(t= seq(0, 15, by = 0.01), phi = 0.3){
  ratio =( 1-pchisq(t, df = 1))/(1 - pf(t/phi, df1= 1, df2= 20))
  ratio
  plot(t, ratio)
  abline(h = 1, col = 2)
}
pvalue( phi =.7)
hist(rchisq(10000, df = 20)/20)

plot(pvalue)

load("U:/R/RA/Data/RFI-newdata/resultsimulation/Model7.Line.Concb.RINa.lneut.llymp.lmono.lbaso.Block/Model7_fit.RData")
dev_sim <- fit$phi.hat.pearson
hist(dev_sim, nclass = 100)
summary(dev_sim)
NB_sim <- fit$NB.disp
hist(NB_sim, nclass = 100)
summary(dev_sim)
plot(y= dev_sim, x= NB_sim)




load("sim_output.RData")
str(sim_output)
