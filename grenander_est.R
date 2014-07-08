library(fdrtool)
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result3/Model102.LineDiet/Model102_result.RData")
str(result$P.values[[3]])
hist(result$P.values[[3]][,"LineDiet"], nclass = 100)

?grenander

z <- result$P.values[[3]][,"LineDiet"]
e <- ecdf(z)
g <- grenander(e)
g
plot(g, log = "y")
plot(g)
plot(density(z))
sum( g$f.knots[-length(g$f.knots)]*diff(g$x.knots) )
install.packages('ADGofTest')
library('ADGofTest')
?ad.test
ad.test(g$F.knots, punif)
plot(g$x.knots, g$F.knots)
sqrt(sum((g$F.knots - g$x.knots)^2))

?ks.test
g$f.knots
ks.test(z, punif)


load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/result4/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")
str(result$P.values[[3]])
ks.test(result$P.values[[3]][,"Blockorder"], punif)
