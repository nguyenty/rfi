runnrep <-function(simdat, criteria){
  nrep <- simdat$nrep
  pi0 <- simdat$pi0
  EE <- simdat$EE
  DE <- simdat$DE
  test.mat.model <- list()
  rt.model <- NULL
  fdr.model <- NULL
  pauc.model <- NULL
  st.model <- NULL
  listout <-  data.frame(Date=as.Date(character()),
                               File=character(), 
                               User=character(), 
                               stringsAsFactors=FALSE) 
    model_th <- 1
    full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                                 neut + lymp + mono + eosi + baso + 
                                 Block + Blockorder)
    #colnames(full_model)
    repeat{
      
      out_model <- fit_model0(full_model, model_th, criteria, simdat)
      
      #vt.model[model_th] <- out_model$vt
      rt.model[model_th] <- out_model$rt
      fdr.model[model_th] <- out_model$fdr
      st.model[model_th] <- out_model$st
      pauc.model[model_th] <- out_model$pauc
      test.mat.model[[model_th]] <- list_model(full_model)$test.mat
      
      assign(paste("ms_criteria", model_th, sep = "_" ),out_model$res_sel)
      proc.time() -pm1
      ms_val <- get(paste("ms_criteria", model_th, sep = "_" ))
      cov_del <- ms_val[1,criteria] # cov_del <- 14; i <- 1
      
      cov_set <- list_model(full_model)$test.mat # dim(cov_set)
      res1 <- data.frame(criteria = colnames(ms_val)[criteria], 
                        model = model_th, 
                        cov_del = rownames(cov_set)[ cov_del])
      listout <- rbind(listout, res1)
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
  ##line effect model######
  model_th <- model_th+1
  full_model <- model.matrix(~Line)
  #colnames(full_model)
  out_model <- fit_model0(full_model, model_th, criteria, simdat)
  rt.model[model_th] <- out_model$rt
  fdr.model[model_th] <- out_model$fdr
  st.model[model_th] <- out_model$st
  pauc.model[model_th] <- out_model$pauc
  test.mat.model[[model_th]] <- list_model(full_model)$test.mat
  
  # oracle model ### 
  
  model_th <- model_th+1
  full_model <- model.matrix(~Line + Concb + RINa + 
                               neut + lymp + mono + 
                               baso + Block)
  #colnames(full_model)
  out_model <- fit_model0(full_model, model_th, criteria, simdat)
  rt.model[model_th] <- out_model$rt
  fdr.model[model_th] <- out_model$fdr
  st.model[model_th] <- out_model$st
  pauc.model[model_th] <- out_model$pauc
  test.mat.model[[model_th]] <- list_model(full_model)$test.mat
#   best.model.est[[nrep]] <- test.mat.model[[which.max(rt.model)]]
#   fdr.est[[nrep]] <- fdr.model
#   rt.est[[nrep]] <- rt.model
   res <- list(best.model.est = 
                 test.mat.model[[which.max(rt.model[-c(model_th, model_th-1)])]], 
               fdr.est = fdr.model,
               rt.est = rt.model, 
               pauc.est = pauc.model, 
               st.est = st.model)
  pi0_dir <- paste0(colnames(ms_val)[criteria], "/pi0_", pi0) #ms_val <- data.frame(pvalue05 = "pvalue05", ks = "ks", ad= "ad", cvm = "cvm")
  dir.create(pi0_dir, showWarnings = FALSE)
   save(res, file = paste0(pi0_dir, "/res_", nrep, ".RData")) # res <-1; 
   write.csv(listout, file = paste0(pi0_dir, "/listout_",nrep,  ".csv"), row.names = FALSE)
}
