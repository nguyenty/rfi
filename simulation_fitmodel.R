fit_model0 <- function(full_model, model_th, criteria, simdat){ # model_th <- 1
  y <- simdat$y
  pi0 <- simdat$pi0
  DE <- simdat$DE
  EE <- simdat$EE
  log.offset <- log(apply(y, 2, quantile, 0.75))
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit2 <- QL.fit(y, design.list, test.mat, # dim(counts)
                 log.offset = log.offset, print.progress=FALSE,
                 Model = "NegBin", method = "optim")
  result2<- QL.results(fit2, Plot = FALSE)
  
  ql <- result2$Q.values[[3]][,"Line"]
  rt <- sum(ql<=.05) # mean(s%in%degene) mean(result2$Q.values[[3]][,"Line"]<.05)
  vt <- sum(ql[1:EE]<=.05) # mean(s%in%degene) mean(result2$Q.values[[3]][,"Line"]<.05)
  fdr <- vt/max(rt, 1) # rt <- 0
  lab <- c(rep(0, EE), rep(1, DE))
  pauc <- pauc_out(ql, lab)
  res_sel <- sel_criteria(result2)  
  print(paste("Model", model_th, sep = " "))
  return(list(res_sel = res_sel, fdr = fdr, rt = rt, st = rt - vt, pauc = pauc ))
}
