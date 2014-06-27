## Covariate Set
covset <- read.csv("covset.csv")
attach(covset)

## Load the result of model 11
load("Model11_result.RData") # use str(result) to see what in it

## Load count data 
scount <- read.table("single end uniquely mapped reads count table for Yet.txt", 
                     header = T)

## List of Genes used to find DE Genes

counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>1 & 
                             rowSums(scount[,-1][,Line ==1] > 0) >0 &
                             rowSums(scount[,-1][, Line ==2] >0) >0 ,-1])

## List of Line DE Genes when FDR is controled at 0.05, 0.10, 0.15

degene05 <- which(result$Q.values[[3]][,"Line"]<=0.05)
degene10 <- which(result$Q.values[[3]][,"Line"]<=0.10)
degene15 <- which(result$Q.values[[3]][,"Line"]<=0.15)

## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15

degene <- c(length(degene05), length(degene10), length(degene15))

## List of subset of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15

lf105 <- which(abs(log(apply(counts[degene05, Line ==1]+1, 1, mean)/apply(counts[degene05, Line ==2]+1, 1, mean))) >=1)

lf110 <- which(abs(log(apply(counts[degene10, Line ==1]+1, 1, mean)/apply(counts[degene10, Line ==2]+1, 1, mean))) >=1)

lf115 <- which(abs(log(apply(counts[degene15, Line ==1]+1, 1, mean)/apply(counts[degene15, Line ==2]+1, 1, mean))) >=1)

## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15

lf1 <- c(length(lf105), length(lf110), length(lf115))

## Summary Table
out <- data.frame(FDR =c(0.05, 0.10, 0.15), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "log(FC)>=1")
xtable(out)
