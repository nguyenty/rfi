source("simulation_loadModel7.R")
source("simulation_criteria.R")
source("simulation_listmodel.R")
source("simulation_fitmodel.R")
source("simulation_simcount.R")
source("simulation_runnrep.R")
## simulation replication #####
pm1 <- proc.time()
size <- 5000
pi0 <- .7
  vnrep <- 13:20
# vnrep <- 34:40
# vnrep <- 57:60
# vnrep <- 75:80
# vnrep <- 96:100

nc <- c("pvalue05", "ad")
criteria <- 2 #pvalue05
for(nrep in vnrep) # vnrep <- 1:2
{ pi0_dir <- paste0(nc[criteria], "/pi0_", pi0) #ms_val <- data.frame(pvalue05 = "pvalue05", ks = "ks", ad= "ad", cvm = "cvm")
  pi0_dir1 <- paste0(nc[1], "/pi0_", pi0) #ms_val <- data.frame(pvalue05 = "pvalue05", ks = "ks", ad= "ad", cvm = "cvm")
  dir.create(pi0_dir, showWarnings = FALSE)
  
  if (criteria ==1){
    simdat <- simcount(nrep, size, pi0)
    save(simdat, file = paste0(pi0_dir,"/simdat_", nrep, ".RData"))  
  }
  if (criteria ==2){
    
    load(file= paste0(pi0_dir1,"/simdat_", nrep, ".RData"))
  }
  print(paste0("nrep = ", nrep ))
  runnrep(simdat, criteria)
}


proc.time()-pm1
# 
# nrep <- 1

# 
# nrep <- 1
# 
# simdat <- simcount(nrep, size)
# sim_output <- simdat
# y <- sim_output$y
# s <- sim_output$s
# log.offset <- log(apply(y, 2, quantile, 0.75))
# 
# model_th <- 1
# full_model <- model.matrix(~Line + Concb + RINa + 
#                              neut + lymp + mono + baso + 
#                              Block)
# 
# list_out <- list_model(full_model)
# design.list <- list_out$design.list
# test.mat <- list_out$test.mat
# fit2 <- QL.fit(y, design.list, test.mat, # dim(counts)
#                log.offset = log.offset, print.progress=FALSE,
#                Model = "NegBin", method = "optim")
# result2<- QL.results(fit2, Plot = FALSE)
# 
# 
# out2 <- table(s%in%degene, result2$Q.values[[3]][,"Line"]<=.05) # mean(s%in%degene) mean(result2$Q.values[[3]][,"Line"]<.05)
# if(dim(out2)[2]==2){
#   rt <- out2[1,2] + out2[2,2]
#   vt <- out2[1,2]
#   fdr <- vt/rt  
# } else{  rt <- 0; fdr <- 0
#          
# }
# out2
# fdr
# rt
# str(result2)
# str(fit2)
# hist(fit2$phi.hat.dev, nclass = 100)
# hist(result2$P.values[[3]][!(s%in%degene), "Line"], nclass = 100)
# hist(result2$P.values[[3]][(s%in%degene), "Line"], nclass = 100)
# hist(result2$P.values[[3]][, "Line"], nclass = 100)
# which(result2$Q.values[[3]][, "Line"] <= .05)
# x <- max(result2$P.values[[3]][which(result2$Q.values[[3]][, "Line"] <= .05), "Line"])
# sum(result2$P.values[[3]][!(s%in%degene), "Line"]<=x)
# simdat$s_degene
# str(fit2)
# pl <- 1-pchisq(fit2$LRT[,1], 1)
# sort(pl)[263]
# 263-sum(pl[(s%in%degene)]<= sort(pl)[263])
# sum(pl[!(s%in%degene)]<= sort(pl)[263])
# 
# hist(pl[!(s%in%degene)], nclass = 100)
# hist(pl[(s%in%degene)], nclass = 100)