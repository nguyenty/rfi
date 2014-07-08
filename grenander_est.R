library(fdrtool)

g_cdf <- function(z){
  e <- ecdf(z)
  g <- grenander(e)
  g
}
sel_criteria <- function(result){
  dat <- result$P.values[[3]][,colnames(result$P.values[[3]])]
  # Crames Von Miser statistics
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
return( data.frame(pvalue05 = order(pvalue_05),
                   ad = order(ad),
                   cvm = order(cvm),
                   ks = order(ks)))
}

# model 1
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")

sel_criteria(result)

# model 2
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model2.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block/Model2_result.RData")
sel_criteria(result)

# model 3
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model3.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model3_result.RData")

sel_criteria(result)
# model 4
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model4.Line.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model4_result.RData")

sel_criteria(result)

# model 5
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model5.Line.RFI.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model5_result.RData")

sel_criteria(result)
# model 6

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model6.Line.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model6_result.RData")
sel_criteria(result)

# model 7
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model7.Line.Concb.RINa.lneut.llymp.lmono.lbaso.Block/Model7_result.RData")
sel_criteria(result)

# model 8
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model8.Line.Concb.RINa.lneut.llymp.lmono.Block/Model8_result.RData")
sel_criteria(result)
