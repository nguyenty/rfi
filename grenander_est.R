library(fdrtool)

# model 1
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
  {
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)

# model 2
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model2.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block/Model2_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
{
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)


# model 3
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model3.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model3_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
{
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)


# model 4
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model4.Line.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model4_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
{
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)



# model 5
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model5.Line.RFI.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model5_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
{
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)

# model 6
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model6.Line.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model6_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
{
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)


# model 7
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model7.Line.Concb.RINa.lneut.llymp.lmono.lbaso.Block/Model7_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
{
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)


# model 8
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model8.Line.Concb.RINa.lneut.llymp.lmono.Block/Model8_result.RData")
# Crames Von Miser distances
cvm <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])], 2, function(z)
{
  e <- ecdf(z)
  g <- grenander(e)
  cvm <- sum((g$F.knots - g$x.knots)^2 *diff(c(0,g$x.knots)))
  cvm
})

which.min(cvm)


# Proportion of pvalue less than 0.05
pvalue_05 <- apply(result$P.values[[3]][,colnames(result$P.values[[3]])]<=0.05, 2, sum)
which.min(pvalue_05)
order(pvalue_05)
order(cvm)

