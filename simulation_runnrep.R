runnrep <-function(simdat, criteria){
  nrep <- simdat$nrep
  test.mat.model <- list()
  rt.model <- NULL
  fdr.model <- NULL
  list_cov_out1 <-  data.frame(Date=as.Date(character()),
                               File=character(), 
                               User=character(), 
                               stringsAsFactors=FALSE) 
    model_th <- 1
    full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                                 neut + lymp + mono + eosi + baso + 
                                 Block + Blockorder)
    #colnames(full_model)
    repeat{
      
      out_model <- fit_model(full_model, model_th, criteria, simdat)
      
      #vt.model[model_th] <- out_model$vt
      rt.model[model_th] <- out_model$rt
      fdr.model[model_th] <- out_model$fdr
      test.mat.model[[model_th]] <- list_model(full_model)$test.mat
      
      assign(paste("ms_criteria", model_th, sep = "_" ),out_model$res_sel)
      proc.time() -pm1
      ms_val <- get(paste("ms_criteria", model_th, sep = "_" ))
      cov_del <- ms_val[1,criteria] # cov_del <- 14; i <- 1
      
      cov_set <- list_model(full_model)$test.mat # dim(cov_set)
      res1 <- data.frame(criteria = colnames(ms_val)[criteria], 
                        model = model_th, 
                        cov_del = rownames(cov_set)[ cov_del])
      list_cov_out1 <- rbind(list_cov_out1, res1)
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
#   best.model.est[[nrep]] <- test.mat.model[[which.max(rt.model)]]
#   fdr.est[[nrep]] <- fdr.model
#   rt.est[[nrep]] <- rt.model
   res <- list(best.model.est = test.mat.model[[which.max(rt.model)]], 
                     fdr.est = fdr.model,
                     rt.est = rt.model)
   save(res, file = paste0(colnames(ms_val)[criteria], "/res_", nrep, ".RData"))
   write.csv(list_cov_out1, file = paste0(colnames(ms_val)[criteria], "/list_cov_out1new_",nrep,  ".csv"), row.names = FALSE)
}
