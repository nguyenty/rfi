
###g_cdf#####
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
