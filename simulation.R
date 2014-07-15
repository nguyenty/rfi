load("Model7_result.RData")
str(result)
load("Model7_fit.RData")
str(fit)
result$m0["QLSpline", ]/12222
?pairs
boxplot(RINa~Line)
hist(fit$coef[,2], nclass = 100)

which(result$Q.values[[3]][,"Line"] <.05)
sum(which(abs(fit$coef[,2])>0.5) %in%which(result$Q.values[[3]][,"Line"] <=0.05))

covset <- read.table("covset.txt")
attach(covset)
table(covset$Line, covset$Block)
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
pairs(USJudgeRatings[1:5], panel = panel.smooth,
      cex = 1.5, pch = 24, bg = "light blue",
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(USJudgeRatings, lower.panel = panel.smooth, upper.panel = panel.cor)

pairs(log(result$P.values[[3]][]), diag.panel = panel.hist,
      upper.panel = panel.cor)

pairs(fit$coef[], diag.panel = panel.hist,
      upper.panel = panel.cor)


covset <- read.table("covset.txt")
str(covset)
str(fit$coef)
colnames(fit$coef)
pairs(fit$coef)
str(fit$LRT)
colnames(fit$LRT)

str(result$P.values[[3]])
pairs(result$P.values[[3]])
# intercept
hist(fit$coef[,1], nclass = 100)

# Line2

hist(fit$coef[,2], nclass = 100)
shapiro.test(fit$coef[,2])

# Concb
hist(fit$coef[,3], nclass = 100)

# RINa
hist(fit$coef[,4], nclass = 100)

# lneut
hist(fit$coef[,5], nclass = 100)
# llymp
hist(fit$coef[,6], nclass = 100)

# lmono
hist(fit$coef[,7], nclass = 100)

# lbaso
hist(fit$coef[,8], nclass = 100)

# Block2 3 4
hist(fit$coef[,9], nclass = 100)
hist(fit$coef[,10], nclass = 100)
hist(fit$coef[,11], nclass = 100)


load("Model5_result.RData")
str(result)
load("Model5_fit.RData")
str(fit)
str(fit$coef)
colnames(fit$coef)

str(fit$LRT)
colnames(fit$LRT)
# intercept
hist(fit$coef[,1], nclass = 100)
# Line2

hist(fit$coef[,2], nclass = 100)
shapiro.test(fit$coef[,2])

# Concb
hist(fit$coef[,3], nclass = 100)

# RFI
hist(fit$coef[,4], nclass = 100)

# lneut
hist(fit$coef[,5], nclass = 100)
# llymp
hist(fit$coef[,6], nclass = 100)

# lmono
hist(fit$coef[,7], nclass = 100)

# lbaso
hist(fit$coef[,8], nclass = 100)

# Block2 3 4
hist(fit$coef[,9], nclass = 100)
hist(fit$coef[,10], nclass = 100)
hist(fit$coef[,11], nclass = 100)
hist(fit$coef[,12], nclass = 100)
hist(fit$coef[,13], nclass = 100)
hist(fit$coef[,14], nclass = 100)
hist(fit$coef[,15], nclass = 100)
hist(fit$coef[,16], nclass = 100)
hist(fit$coef[,17], nclass = 100)
covset <- read.table("covset.txt")
str(covset)
str(fit$NB.disp)
hist(fit$NB.disp, nclass = 50)
