load("Model7_result.RData")
load("Model7_fit.RData")
result$m0["QLSpline", ]/12222

scount <- read.table("single end uniquely mapped reads count table for Yet.txt", 
                     header = T)

## List of Genes used to find DE Genes

counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>8 ,-1])

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
# check if the condition of count avarage is sastified
# change to simulation setup
#estimate the number 
# the following function is the same as function pval.hist
# in the paper of Pounds et al. 2012, EBT, with the estimation of 
# density obtained from the paper by Korbinian Strimmer 2008
# 
library(fdrtool)
pval.hist.grenander <- function(p.value){
  grenander.out <- grenander(ecdf(p.value))
  p.brks <- c(0, grenander.out$x.knots)
  b.edf <- c(0, grenander.out$F.knots)
  p.diffs <- diff(p.brks)
  h.cdf <- approx(p.brks, b.edf, xout = p.value)$y  # get the histogram-based CDF estimates from each p-value
  p.hist <- exp(log(diff(b.edf))-log(diff(p.brks))) # get the hight for each histogram bar
  pi0.hat <- min(p.hist)                            # get the pi0 estimate from histogram bar
  h.ebp <- approx(p.brks, pi0.hat/c(p.hist, p.hist[length(p.hist)]), xout = p.value)$y # get the conservative EBP interpolation 
  h.fdr <- exp(log(pi0.hat) + log(p.value) - log(h.cdf))                                     # Get the histogram based FDR estimate
  h.ebp[p.value==0] <- 0
  h.fdr[p.value==0] <- 0
  return(list( p.value = p.value,          # input p-value,
               h.cdf = h.cdf,              # the histogram Grenander based cdf estimate
               h.fdr = h.fdr,              # the histogram Grenander based FDR
               h.ebp = h.ebp,              # the histogram Grenander based EBP
               p.brks = p.brks,            # the p-value break-points of the histogram
               p.hist = p.hist,            # the heights of each histogram bar
               edf.brks = b.edf,           # the breaks points in the EDF of the histogram estimator
               pi0.hat = pi0.hat))         # the histogram Grenander based estimate of the proportion of tests with a true null hypothesis
}

pvalue_line <- result$P.values[[3]][,"Line"]
qvalue_line <- result$Q.values[[3]][,"Line"]
gre_out <- pval.hist.grenander(pvalue_line)
ebp_line <- gre_out$h.ebp
hist(ebp_line, nclass = 100)
hist(pvalue_line, nclass = 100)
length(pvalue_line)
result$m0["QLSpline", "Line"]
sum(ebp_line >0.5)/length(pvalue_line)
plot(ebp_line, pvalue_line)


hist(pvalue_line[ebp_line<0.5], nclass = 100)
hist(fit$coef[,2][ebp_line<0.5], nclass = 100)
hist(fit$coef[,2], nclass = 100)
summary(fit$coef[,2])
boxplot(fit$coef[,2])
sum(fit$coef[,2][ebp_line<0.5] >0.06)
sum(fit$coef[,2] >0.06)
coef_line<- fit$coef[,2]
plot(ebp_line, coef_line)

fdr_line <- gre_out$h.fdr
sum(fdr_line <=.10)
which ()
plot(pvalue_line, fdr_line)
plot(pvalue_line, qvalue_line)
plot(fdr_line, qvalue_line, type = "l")
s <- sample(length(pvalue_line), 10000)
ind <- s[order(s)]
spvalue_line <- pvalue_line[s]
sebp_line <- pval.hist.grenander(spvalue_line)$h.ebp
sum(sebp_line>0.5)/length(spvalue_line)
full_model <- model.matrix(~Line + Concb + RINa + lneut + llymp + lmono + lbaso + Block)
Block <- as.factor(Block)
Line <- as.factor(Line)
coef_beta <- fit$coef 
coef_beta[,2] <- fit$coef[,2]*(ebp_line<0.5)
dim(coef_beta)
colnames(full_model)
set.seed(1)
J <- 5000
s <- sample(dim(coef_beta)[1], J)
s <- s[order(s)]
mu <- coef_beta%*%t(full_model)
omega <- fit$NB.disp
degene <- which((ebp_line<0.5))
s_mu <- mu[s,]
s_omega <- omega[s]

degene <- which(ebp_line<0.5)

s_degene <- s[s%in%degene]  

