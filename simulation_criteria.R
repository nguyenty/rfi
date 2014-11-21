
###g_cdf#####
g_cdf <- function(z){
  e <- ecdf(z)
  g <- grenander(e)
  g
}
sel_criteria <- function(result){
  dat <- result$P.values[[3]][,colnames(result$P.values[[3]])]
  dat2 <- result$Q.values[[3]][,colnames(result$Q.values[[3]])]
  # Crames Von Miser statistics
  if(is.vector(dat)) dat <- as.matrix(dat)
  if(is.vector(dat2)) dat2 <- as.matrix(dat2)
#   cvm <- apply(dat, 2, function(z)sum((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2 *
#                                         diff(c(0,g_cdf(z)$x.knots))))
#   # Kolmogorow Smirnov statistics 
#   ks <- apply(dat, 2, function(z)max((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2))
#   
  # Anderson-Darling statistics
  ad <- apply(dat, 2, function(z)sum((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2/
                                       g_cdf(z)$x.knots*(1-g_cdf(z)$x.knots)*
                                       diff(c(0,g_cdf(z)$x.knots))))
  # Proportion of pvalue less than 0.05
  pvalue_05 <- apply(dat<=0.05, 2, sum)
#   qvalue05 <- apply(dat2<=0.05, 2, sum)
  out <- data.frame(pvalue05 = order(pvalue_05),
                    ad = order(ad))
  return(out)
}

pauc_out <- function(test.vector, lab){
  lab <- as.factor(lab)
  roc.out <- roc(1-test.vector, lab) # plot(roc.out)
  roc.ind <- sum(roc.out$fpr<=.05)
  roc.min <- roc.out$cutoffs[roc.ind]
  pauc <- auc(roc.out, min =roc.min)
  return(pauc)
}
