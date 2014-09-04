load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsingle/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")
#str(result)
Model1 <- result$P.values[[3]][, "Line"]
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsingle/Model2.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block/Model2_result.RData")
Model2 <- result$P.values[[3]][, "Line"]


load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsingle/Model3.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model3_result.RData")
Model3 <- result$P.values[[3]][, "Line"]


load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsingle/Model4.Line.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model4_result.RData")

Model4 <- result$P.values[[3]][, "Line"]

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsingle/Model5.Line.RFI.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model5_result.RData")
Model5 <- result$P.values[[3]][, "Line"]

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsingle/Model6.Line.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model6_result.RData")
Model6 <- result$P.values[[3]][, "Line"]

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultsingle/Model7.Line.Concb.RINa.lneut.llymp.lmono.lbaso.Block/Model7_result.RData")
Model7 <- result$P.values[[3]][, "Line"]

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


dat <- cbind(Model1, Model2,
             Model3,Model4,
             Model5,Model6,
             Model7)
colnames(dat) <- paste("M", 1:7, sep = "")
pairs(dat, upper.panel = panel.cor,
      cex = 1.5, pch = 24, bg = "light blue",
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2,
      main = "Single End Data")



#################paired End data##################

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")
#str(result)
Model1 <- result$P.values[[3]][, "Line"]
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired/Model2.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block/Model2_result.RData")
Model2 <- result$P.values[[3]][, "Line"]


load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired/Model3.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model3_result.RData")
Model3 <- result$P.values[[3]][, "Line"]


load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired/Model4.Line.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.lbaso.Block/Model4_result.RData")

Model4 <- result$P.values[[3]][, "Line"]

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired/Model5.Line.RFI.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model5_result.RData")
Model5 <- result$P.values[[3]][, "Line"]

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired/Model6.Line.Concb.RINb.RINa.lneut.llymp.lmono.lbaso.Block/Model6_result.RData")
Model6 <- result$P.values[[3]][, "Line"]

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/RFI-newdata/resultpaired/Model7.Line.Concb.RINa.lneut.llymp.lmono.lbaso.Block/Model7_result.RData")
Model7 <- result$P.values[[3]][, "Line"]

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


dat <- cbind(Model1, Model2,
             Model3,Model4,
             Model5,Model6,
             Model7)
colnames(dat) <- paste("M", 1:7, sep = "")
pairs(dat, upper.panel = panel.cor,
      cex = 1.5, pch = 24, bg = "light blue",
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2,
      main = "Paired End Data")