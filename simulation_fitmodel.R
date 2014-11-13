fit_model <- function(full_model, model_th, criteria, sim_output){ # model_th <- 1
  y <- sim_output$y
  s <- sim_output$s
  log.offset <- log(apply(y, 2, quantile, 0.75))
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit2 <- QL.fit(y, design.list, test.mat, # dim(counts)
                 log.offset = log.offset, print.progress=FALSE,
                 Model = "NegBin", method = "optim")
  result2<- QL.results(fit2, Plot = FALSE)
  
  
  out2 <- table(s%in%degene, result2$Q.values[[3]][,"Line"]<=.05) # mean(s%in%degene) mean(result2$Q.values[[3]][,"Line"]<.05)
  if(dim(out2)[2]==2){
    rt <- out2[1,2] + out2[2,2]
    vt <- out2[1,2]
    fdr <- vt/rt  
  } else{  rt <- 0; fdr <- 0
           
  }
  res_sel <- sel_criteria(result2)  
  print(paste("Model", model_th, sep = " "))
  return(list(res_sel = res_sel, fdr = fdr, rt = rt))
}
