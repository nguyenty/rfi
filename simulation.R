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
fdr_line <- gre_out$h.fdr
Block <- as.factor(Block)
Line <- as.factor(Line)

full_model <- model.matrix(~Line + Concb + RINa + lneut + llymp + lmono + lbaso + Block)
dim(full_model)
coef_beta <- fit$coef 
coef_beta[,2] <- fit$coef[,2]*(ebp_line<0.5)
set.seed(1)
J <- 1000
s <- sample(dim(coef_beta)[1], J)
s <- s[order(s)]

offset <- apply(counts, 2, quantile, 0.75)
lf=length(offset)

R=matrix(rep(offset,dim(coef_beta)[1]),
         ncol=lf,
         byrow=T)


mu <- exp(coef_beta%*%t(full_model))* R
omega <- fit$NB.disp
degene <- which((ebp_line<0.5))

s_mu <- mu[s,]
s_omega <- omega[s]
s_degene <- intersect(degene, s) 


##Sim counts data
y <- array(0, dim = c(J,31))

for(j in 1:J){
  repeat{
      for(k in 1:31){
        y[j,k]  <- rnbinom(n=1, size=1/s_omega[j], mu=s_mu[j,k])
      }
      if (mean(y[j,])>8& sum(y[j,]>0)>3) break
    }
}


