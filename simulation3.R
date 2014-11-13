library(Matrix)
source("simulation_loadModel72.R")
source("simulation_criteria.R")
source("simulation_listmodel.R")

## Function do all the things with input Full model
# sim_output <- list(used_gene = used_gene, 
#                    used_omega = used_omega, 
#                    used_count = used_count, 
#                    degene_original = degene, 
#                    s=s, s_degene = s_degene, 
#                    y = y)
fit_model <- function(full_model, model_th, criteria, sim_output){ # model_th <- 1
  y <- sim_output$y
  s <- sim_output$s
  log.offset <- log(apply(y, 2, quantile, 0.75))
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit2 <- QL.fit(y, design.list, test.mat, # dim(counts)
                log.offset = log.offset, print.progress=FALSE,
                Model = "NegBin")
  result2<- QL.results(fit2, Plot = FALSE)
  
  fit <- QL.fit(y, design.list, test.mat, # dim(counts)
                 log.offset = log.offset, print.progress=FALSE,
                 Model = "NegBin")
  result<- QL.results(fit, Plot = FALSE)
  
  
  
  #result2 <- result out2 <- out
  out <- table(s%in%degene, result$Q.values[[3]][,"Line"]<=.05) # mean(s%in%degene) mean(result2$Q.values[[3]][,"Line"]<.05)
  out2 <- table(s%in%degene, result2$Q.values[[3]][,"Line"]<=.05) # mean(s%in%degene) mean(result2$Q.values[[3]][,"Line"]<.05)
  if(dim(out2)[2]==2){
    rt <- out2[1,2] + out2[2,2]
    vt <- out2[1,2]
    fdr <- vt/rt  
  } else{  rt <- 0; fdr <- 0
    
  }
  res_sel <- sel_criteria(result2)  
  k <- nrow(test.mat)
  name_model <- NULL 
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(resultdir,"/", colnames(res_sel)[criteria], "/Model",model_th,name_model, sep ="")
  dir.create(model_dir, showWarnings = FALSE)
  save(result2, file = paste(model_dir,"/Model",model_th, "_result.RData", sep =""))
  save(fit2, file = paste(model_dir,"/Model",model_th, "_fit.RData", sep =""))
  for(i in 1:(nrow(test.mat))){
    postscript(paste(model_dir,"/Model", 
                     model_th, row.names(test.mat)[i],".eps", sep =""))
    hist(result2$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
    
    pdf(paste(model_dir,"/Model", 
              model_th, row.names(test.mat)[i],".pdf", sep =""))
    hist(result2$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
  }
  print(paste("Model", model_th, sep = " "))
  return(list(res_sel = res_sel, fdr = fdr, rt = rt))
}


fdr.est <- rt.est <- best.model.est <- list()

## simulation replication #####
pm1 <- proc.time()
for(nrep in c(1)) # nrep <- 500
{
  test.mat.model <- list()
  #vt.model <- NULL
  rt.model <- NULL
  fdr.model <- NULL
J <- dim(used_beta)[1]
#used_gene is the list of indexes of genes from the original data 12280 
size <- 8000
set.seed(nrep)
s <- sample( dim(used_beta)[1], size = size)
s <- s[order(s)]
s_mu <- mu[s,]
s_omega <- used_omega[s]
s_degene <- intersect(degene, s) 
#length(s_degene)/length(s)

##Sim counts data

y <- array(0, dim = c(size,31))

# nrep <- 100;size <- 5000 ; j <-5000; k <- 1

for(j in 1:size){ # j <- 1; k <- 1
  repeat{
    for(k in 1:31){
      #set.seed((32*nrep+k)*size + j)
      y[j,k]  <- rnbinom(n=1, size=1/s_omega[j], mu=s_mu[j,k])
    }
    if (mean(y[j,])>8& sum(y[j,]>0)>3) break
  }
}



# nrep <- 1
sim_outputnew <- list(used_gene = used_gene, used_omega = used_omega, 
                   used_count = used_count, 
                   degene_original = degene, 
                   s = s,  s_degene = s_degene, 
                   y = y, nrep = nrep)
save(sim_outputnew, file = paste0("sim_outputnew_", nrep, ".RData"))



list_cov_out1 <-  data.frame(Date=as.Date(character()),
                                 File=character(), 
                                 User=character(), 
                                 stringsAsFactors=FALSE) 

for(criteria in 1){ # i <- 1
  model_th <- 1
  full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                               neut + lymp + mono + eosi + baso + 
                               Block + Blockorder)
  #colnames(full_model)
  repeat{
    
    out_model <- fit_model(full_model, model_th, criteria, sim_outputnew)
    
    #vt.model[model_th] <- out_model$vt
    rt.model[model_th] <- out_model$rt
    fdr.model[model_th] <- out_model$fdr
    test.mat.model[[model_th]] <- list_model(full_model)$test.mat
    
    assign(paste("ms_criteria", model_th, sep = "_" ),out_model$res_sel)
    proc.time() -pm1
    ms_val <- get(paste("ms_criteria", model_th, sep = "_" ))
    cov_del <- ms_val[1,criteria] # cov_del <- 14; i <- 1
    
    cov_set <- list_model(full_model)$test.mat # dim(cov_set)
    res <- data.frame(criteria = colnames(ms_val)[criteria], 
                      model = model_th, 
                      cov_del = rownames(cov_set)[ cov_del])
    list_cov_out1 <- rbind(list_cov_out1, res)
    if (cov_del ==1) break
    block_ind <- grep("Block2", colnames(full_model))
    blockorder_ind <-grep("Blockorder", colnames(full_model))
    indicator <- FALSE
    if(length(block_ind)!=0){
      if(cov_del+1 == cov_set["Block",2])
      {
        full_model <- full_model[, -c(block_ind, block_ind+1, block_ind+2)]
        indicator <- TRUE
      }
    }
    
    if(length(blockorder_ind)!=0){
      if(cov_del+1 == cov_set["Blockorder", 2])
      {
        full_model <- full_model[, -blockorder_ind] 
        indicator <- TRUE
      }
    }
    
    if (indicator == FALSE) full_model <- full_model[, -(cov_del +1)] 
    model_th <- model_th +1
  }
}
 
  best.model.est[[nrep]] <- test.mat.model[[which.max(rt.model)]]
  fdr.est[[nrep]] <- fdr.model
  rt.est[[nrep]] <- rt.model
res.outnew <- list(best.model.est = best.model.est[[nrep]], 
                fdr.est = fdr.est[[nrep]],
                rt.est = rt.est[[nrep ]])
# save(res.outnew, file = paste0("res.outnew_", nrep, ".RData"))
#  write.csv(list_cov_out1, file = paste0("list_cov_out1new_",nrep,  ".csv"), row.names = FALSE)
 #write.csv(list_cov_out1, file = paste0("list_cov_out2_",nrep,  ".csv"), row.names = TRUE)
#  write.csv(best.model, file = paste0("best_model1_", nrep, ".csv"), row.names = TRUE)
#  write.csv(best.model, file = paste0("best_model2_", nrep, ".csv"), row.names = TRUE)
}

proc.time()-pm1
# 
# fdr.est
# rt.est
# best.model
# 
# ## combine all results together ######
# nrep <- 44
# res.all <- list()
# 
# nrep
# fdr.sim <- list()
# rt.sim <- list()
# best.model.sim <- NULL
# for( i in 1:nrep){
#   path <- paste0("res.outnew_", i, ".RData")
#   load(path)
#   res.all[[i]] <- res.outnew
#   best.model.sim[i] <- which.max(res.outnew$rt)
#   rt.sim[[i]] <- res.outnew$rt
#   fdr.sim[[i]] <- res.outnew$fdr
# }
# 
# rt.sim[[1]]
# best.model.sim
# str(res.all[[8]])
# rt.best.sim <- NULL
# fdr.best.sim <- NULL
# for(i in 1:nrep){
#   fdr.best.sim[i] <- fdr.sim[[i]][which.max(rt.sim[[i]])]
#   rt.best.sim[i] <- max(rt.sim[[i]])
# }
# fdr.best.sim
# fdr.mean <- mean(fdr.best.sim)
# fdr.mean
# fdr.se <- sd(fdr.best.sim)/sqrt(nrep)
# fdr.se
# rt.best.sim
# 
# 
# ### run the model with only Line effect#####
# st.line <-rt.line <- fdr.line <- NULL
# for(i in 1:44){ # i <- 1
#   path <- paste0("sim_outputnew_", i, ".RData")
#   load(path)
#   model_th <- 100+i# i <- 1
#   full_model <- model.matrix(~Line)
#   out_model <- fit_model(full_model, model_th, 1, sim_outputnew)
#   rt.line[i] <- out_model$rt
#   fdr.line[i] <- out_model$fdr
#   st.line[i] <- rt.line[i] - round(rt.line[i]*fdr.line[i])
# }
# sim.line <- list(rt.line = rt.line, fdr.line = fdr.line, st.line = st.line)
# save(sim.line, "sim.line.RData")
# rt.line
# fdr.line
# st.line
# 
# library(AUC)
# label <- as.factor(as.numeric((s%in%degene)))
# auc(roc(1-result2$Q.values[[3]][, "Line"], label))
# auc(roc(1-result$Q.values[[3]], label))

auc_out <- function(test.vector, lab){
  lab <- as.factor(lab)
  roc.out <- roc(1-test.vector, lab) # plot(roc.out)
  roc.ind <- sum(roc.out$fpr<=.05)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- auc(roc.out, min =roc.min)
  return(pauc)
}

auc_out(result2$Q.values[[3]][, "Line"], label)
auc_out(result$Q.values[[3]], label)
