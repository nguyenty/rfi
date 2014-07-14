load("Model7_result.RData")
str(result)
load("Model7_fit.RData")
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