#install.packages("rjags")
library(rjags)
set.seed(20141020)
#############modelm - using point mass mixture prior for signals###############
modelm <- "
model{
# likelihood 
for (i in 1:length(y)){
y[i] ~ dnorm(beta[i], 1/sigma[i]^2)
beta[i] <- (1 - bin.beta[i])*norm.beta[i]
bin.beta[i] ~ dbern(pi.beta)
norm.beta[i] ~ dnorm(beta0, 1/sigma0^2)
sigma[i] ~ dunif(0, 100)
}

# prior distribution for the parameters ####
pi.beta ~ dbeta(5, 1)
beta0 ~ dnorm(0, 100)
sigma0 ~ dunif(0, 100)
}
"
# load("U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/ks/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_fit.RData")
# load("U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/ks/Model1.Line.Diet.RFI.Concb.RINb.Conca.RINa.lneut.llymp.lmono.leosi.lbaso.Block.Blockorder/Model1_result.RData")

load("Model1_fit.RData")
load("Model1_result.RData")

data <- list(y = fit$coef[,3])

m0 <- proc.time()
mm <- jags.model(textConnection(modelm), data,n.chains = 1) # mix point mass
resm <- coda.samples(mm, c("beta","sigma","beta0","sigma0","pi.beta", 
                           "bin.beta")
                           , 5000) # mix point mass


proc.time() -m0
which(colnames(resm[[1]])=="sigma[1]")
which(colnames(resm[[1]])=="sigma[12280]")
dim(resm[[1]])
resm[[1]][, "sigma[1]"]
sigma.i <- apply(resm[[1]][,24563:36842], 2, mean)
sigma.0 <- mean(resm[[1]][, "sigma0"])
beta.0 <- mean(resm[[1]][,"beta0"])
summary(sigma.i)
hist(sigma.i, nclass = 100)
abline(v = sigma.0)
#resm[[1]][,"beta[1]"]
new.beta <- fit$coef[,3]*(1/sigma.i^2)/(1/sigma.i^2 + 1/sigma.0^2) + 
                            beta.0 * (1/sigma.0^2)/(1/sigma.i^2 + 1/sigma.0^2)

hist(fit$coef[,3], nclass = 1000, xlim = c(-5,5))

hist(new.beta, nclass = 1000)
hist(fit$coef[,3], nclass = 1000)
resm[[1]][,"beta[1]"]
resm[[1]][,"bin.beta[1]" ]
save(new.beta, file = "new.beta.RData")
