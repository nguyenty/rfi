counts <- read.csv("counts_12280.csv")
dir <- "U:/R/RA/Data/RFI-newdata/resultpairedlogcbc/pvalue05/"
#dir <- "U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue05/"
#dir <- "U:/R/RA/Data/RFI-newdata/resultsimulation/pvalue05/"
f <- list()
r <- list()
q.line <- vector()
dir.list <- list.files(dir)
for(i in 1:length(dir.list)){ # i <- 7
  diri <- paste0(dir, dir.list[grep(paste0("Model", i,".Line" ), dir.list)], "/")
  load(file = paste0(diri, paste0("Model", i, "_fit.RData")))
  load(file = paste0(diri, paste0("Model", i, "_result.RData")))
  f[[i]] <- fit
  r[[i]] <- result
  q.line[i] <- sum(r[[i]]$Q.values[[3]][,"Line"]<0.05)
  
}

### genes list including gene ID log2FC, pvalue, qvalue of those genes ######

log2fc.line.logcbc <- abs(log2(exp(f[[7]]$coef[,2])))
pvalue.line.logcbc <- r[[7]]$P.values[[3]][, "Line"]
qvalue.line.logcbc <- r[[7]]$Q.values[[3]][, "Line"]



dir <- "U:/R/RA/Data/RFI-newdata/resultpairedcbc/pvalue05/"

f2 <- list()
r2 <- list()
q2.line <- vector()
dir.list <- list.files(dir)
for(i in 1:length(dir.list)){ # i <- 7
  diri <- paste0(dir, dir.list[grep(paste0("Model", i,".Line" ), dir.list)], "/")
  load(file = paste0(diri, paste0("Model", i, "_fit.RData")))
  load(file = paste0(diri, paste0("Model", i, "_result.RData")))
  f2[[i]] <- fit
  r2[[i]] <- result
  q2.line[i] <- sum(r2[[i]]$Q.values[[3]][,"Line"]<0.05)
  
}


log2fc.line.cbc <- abs(log2(exp(f2[[7]]$coef[,2])))
pvalue.line.cbc <- r2[[7]]$P.values[[3]][, "Line"]
qvalue.line.cbc <- r2[[7]]$Q.values[[3]][, "Line"]

GeneList <- as.data.frame(cbind( 
      log2fc.line.logcbc = log2fc.line.logcbc, 
      pvalue.line.logcbc = pvalue.line.logcbc,
      qvalue.line.logcbc = qvalue.line.logcbc,
      log2fc.line.cbc = log2fc.line.cbc, 
      pvalue.line.cbc = pvalue.line.cbc,
      qvalue.line.cbc = qvalue.line.cbc))

rownames(GeneList) <- counts[,1]
write.csv(GeneList, "GeneList.csv")
dim(GeneList[(GeneList$qvalue.line.cbc<=0.05),])

dim(GeneList[(GeneList$log2fc.line.logcbc>=1)&
               (GeneList$qvalue.line.logcbc<=0.1),])

dim(GeneList[(GeneList$log2fc.line.cbc>=1)&
               (GeneList$qvalue.line.cbc<=0.05),])

dim(GeneList[(GeneList$log2fc.line.logcbc>=1)& 
               (GeneList$qvalue.line.logcbc<=0.1)&
               (GeneList$log2fc.line.cbc>=1)& 
               (GeneList$qvalue.line.cbc<=0.1),])

