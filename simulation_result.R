
## combine all results together ######
nrep <- 100
res.all <- list()

nrep
fdr.sim <- list()
rt.sim <- list()
best.model.sim <- NULL
list.out <- list()
for( i in 1:nrep){ # i <- 1
  path <- paste0("pvalue05/res_", i, ".RData")
  path2 <- paste0("pvalue05/list_cov_out1new_", i, ".csv")
  load(path)
  covout <- read.csv(path2)
  list.out[[i]] <-covout$cov_del[1:6]
  res.all[[i]] <- res
  best.model.sim[i] <- which.max(res$rt)
  rt.sim[[i]] <- res$rt
  fdr.sim[[i]] <- res$fdr
  
}
fdr.sim
rt.sim[[1]]
best.model.sim
str(res.all[[8]])
rt.best.sim <- NULL
fdr.best.sim <- NULL
for(i in 1:nrep){
  fdr.best.sim[i] <- fdr.sim[[i]][which.max(rt.sim[[i]])]
  rt.best.sim[i] <- max(rt.sim[[i]])
}
fdr.best.sim
fdr.mean <- mean(fdr.best.sim)
fdr.mean
fdr.se <- sd(fdr.best.sim)/sqrt(nrep)
fdr.se
rt.best.sim
rt.best.mean <- mean(rt.best.sim)
rt.best.se <- sd(rt.best.sim)/sqrt(nrep)
rt.best.mean
rt.best.se
mean(best.model.sim == 7)

ind <- which(best.model.sim == 7)
ind
fdr.m <- NULL
for(i in 1:length(ind))
{fdr.m[i] <- fdr.best.sim[[ind[i]]]}
mean(fdr.m)

list.out[[ind[1]]]
list.out[ind[1]]
d1 <- NULL
d2 <- as.character(list.out[[ind[1]]])
for (i in 1:length(ind)){
  d1 <- union(d1, as.character(list.out[[ind[i]]]))
  d2 <- intersect(d2, as.character(list.out[[ind[i]]]))
}
d1
d2

#####cvm#####

## combine all results together cvm criterion ######
nrep <- 100
res.all <- list()

nrep
fdr.sim <- list()
rt.sim <- list()
list.out <- list()
best.model.sim <- NULL
for( i in 1:nrep){
  path <- paste0("cvm/res_", i, ".RData")
  path2 <- paste0
  load(path)
  res.all[[i]] <- res
  best.model.sim[i] <- which.max(res$rt)
  rt.sim[[i]] <- res$rt
  fdr.sim[[i]] <- res$fdr
}
fdr.sim
rt.sim[[1]]
best.model.sim
str(res.all[[8]])
rt.best.sim <- NULL
fdr.best.sim <- NULL
for(i in 1:nrep){
  fdr.best.sim[i] <- fdr.sim[[i]][which.max(rt.sim[[i]])]
  rt.best.sim[i] <- max(rt.sim[[i]])
}
fdr.best.sim
fdr.mean <- mean(fdr.best.sim)
fdr.mean
fdr.se <- sd(fdr.best.sim)/sqrt(nrep)
fdr.se
rt.best.sim
mean(best.model.sim == 7)
####

### run the model with only Line effect#####
st.line <-rt.line <- fdr.line <- NULL
for(i in 1:nrep){ # i <- 1
  path <- paste0("simdat/simdat_", i, ".RData")
  load(path)
  model_th <- 100+i# i <- 1
  full_model <- model.matrix(~Line)
  out_model <- fit_model(full_model, model_th, 1, simdat)
  rt.line[i] <- out_model$rt
  fdr.line[i] <- out_model$fdr
  st.line[i] <- rt.line[i] - round(rt.line[i]*fdr.line[i])
}
sim.line <- list(rt.line = rt.line, fdr.line = fdr.line, st.line = st.line)
save(sim.line, file = "simdat/sim.line.RData")
rt.line
fdr.line
st.line
mean(fdr.line)
sd(fdr.line)/10
library(AUC)
label <- as.factor(as.numeric((s%in%degene)))
auc(roc(1-result2$Q.values[[3]][, "Line"], label))
auc(roc(1-result$Q.values[[3]], label))



###result from real data cvm####

dir <- "U:/R/RA/Data/RFI-newdata/resultpairedcbc/cvm/"

f2 <- list()
r2 <- list()
q2.line <- vector()
dir.list <- list.files(dir)
for(i in 1:length(dir.list)){ # i <- 7777
  diri <- paste0(dir, dir.list[grep(paste0("Model", i,".Line" ), dir.list)], "/")
  load(file = paste0(diri, paste0("Model", i, "_fit.RData")))
  load(file = paste0(diri, paste0("Model", i, "_result.RData")))
  f2[[i]] <- fit
  r2[[i]] <- result
  q2.line[i] <- sum(r2[[i]]$Q.values[[3]][,"Line"]<=0.05)
  
}
q2.line

## result from pvalue05#####
dir <- "U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue05/"

f2 <- list()
r2 <- list()
q2.line <- vector()
dir.list <- list.files(dir)
for(i in 1:length(dir.list)){ # i <- 7777
  diri <- paste0(dir, dir.list[grep(paste0("Model", i,".Line" ), dir.list)], "/")
  load(file = paste0(diri, paste0("Model", i, "_fit.RData")))
  load(file = paste0(diri, paste0("Model", i, "_result.RData")))
  f2[[i]] <- fit
  r2[[i]] <- result
  q2.line[i] <- sum(r2[[i]]$Q.values[[3]][,"Line"]<=0.05)
  
}
q2.line
path3 <- paste0("pvalue05/sim.line.RData")
path3

load(path3)
m.rt <- mean(sim.line$rt)
se.rt <- sd(sim.line$rt)/sqrt(nrep)
m.fdr <- mean((sim.line$fdr))
se.fdr <- sd(sim.line$fdr)/sqrt(nrep)
m.rt
se.rt
m.fdr
se.fdr


####pi0 the same size 5000 qvalue05 method######

###result from real data cvm####

pi0 <- c(0.9, 0.8, 0.7, 0.6, 0.5)
i <- 3
dir <- paste0("pvalue05/pi0_", pi0[i])

dir.list <- list.files(dir)
indres <- grep("res", dir.list)
best.model.est <- list()
fdr.est <- list()
rt.est <- list()
for(i in 1:length(indres)){ # i <- 1
  diri <- paste0(dir, "/", dir.list[indres[i]])
  load(file = diri)
  best.model.est[[i]] <- res$best.model.est
  fdr.est[[i]] <- res$fdr.est
  rt.est[[i]] <- res$rt.est
  
}



model.best.sim <- NULL
rt.best.sim <- NULL
fdr.best.sim <- NULL
for(i in 1:length(indres)){
  model.best.sim[i] <- which.max(rt.est[[i]])
  fdr.best.sim[i] <- fdr.est[[i]][model.best.sim[i]]
  rt.best.sim[i] <- max(rt.est[[i]])
}
#fdr.best.sim
fdr.mean <- mean(fdr.best.sim)
fdr.mean
fdr.se <- sd(fdr.best.sim)/sqrt(nrep)
fdr.se
rt.best.sim
model.best.sim
mean(model.best.sim == 7)

## result from pvalue05#####
dir <- "U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue05/"

f2 <- list()
r2 <- list()
q2.line <- vector()
dir.list <- list.files(dir)
for(i in 1:length(dir.list)){ # i <- 7777
  diri <- paste0(dir, dir.list[grep(paste0("Model", i,".Line" ), dir.list)], "/")
  load(file = paste0(diri, paste0("Model", i, "_fit.RData")))
  load(file = paste0(diri, paste0("Model", i, "_result.RData")))
  f2[[i]] <- fit
  r2[[i]] <- result
  q2.line[i] <- sum(r2[[i]]$Q.values[[3]][,"Line"]<=0.05)
  
}
q2.line
path3 <- paste0("pvalue05/sim.line.RData")
path3

load(path3)
m.rt <- mean(sim.line$rt)
se.rt <- sd(sim.line$rt)/sqrt(nrep)
m.fdr <- mean((sim.line$fdr))
se.fdr <- sd(sim.line$fdr)/sqrt(nrep)
m.rt
se.rt
m.fdr
se.fdr