source("simulation_loadModel72.R")
source("simulation_criteria.R")
source("simulation_listmodel.R")
source("simulation_fitmodel.R")
source("simulation_simcount.R")
source("simulation_runnrep.R")
## simulation replication #####
pm1 <- proc.time()
size <- 5000
#  vnrep <- 1:10
#  vnrep <- 11:20
# vnrep <- 21:30
# vnrep <- 31:40
# vnrep <- 41:50
# vnrep <- 51:60
# vnrep <- 61:70
# vnrep <- 71:80
#  vnrep <- 81:90
vnrep <- 91:100

criteria <- 1 #pvalue05
for(nrep in vnrep) # nrep <- 500
{
  simdat <- simcount(nrep, size)
  save(simdat, file = paste0("simdat2/simdat_", nrep, ".RData"))
  print(paste0("nrep = ", nrep ))
  runnrep(simdat, criteria)
}


proc.time()-pm1
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