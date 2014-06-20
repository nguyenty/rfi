require(Matrix)
library(QuasiSeq)
library(edgeR)
### Reading data #######
scount <- read.table("/home/ntyet/research/RFI-newdata/Data for Yet/single end uniquely mapped reads count table for Yet.txt", 
                     header = T)
cbc <- read.table('/home/ntyet/research/RFI-newdata/Data for Yet/CBC data for pigs with RNA-seq data avaible.txt',
                  header =T)

metadata <- read.table("/home/ntyet/research/RFI-newdata/Data for Yet/meta_data_RNA-seq_G9P2.txt", 
                       header = T)

rfiadj <- read.csv("/home/ntyet/research/RFI-newdata/Data for Yet/g8p2_g9p2-rfiadj-FINAL_Jan_2014_rfiadjusted.csv", 
                     header = T)
##### cleaning data####

cbc <- cbc[order(cbc$ear), ]
metadata <- metadata[order(metadata$idpig), ]
rfiadj <- rfiadj[order(rfiadj$idpig),]

fullidpig <- as.numeric(paste("20900", metadata$idpig, sep = ""))

covset <- cbind(metadata[, -4], rfiadj[rfiadj$idpig %in% fullidpig, c("rfi.ADJUSTED")], 
                cbc[, c("iddam", "idsire", "Neutrophil", "Lymphocyte", "Monocyte",
                        "Eosinophil", "Basophil" )])
colnames(covset) <- c("idpig", "Line", "Diet",  "Block", "Blockorder", "Concb", 
                      "RINb", "Conca", "RINa", "RFI",
                      "iddam", "idsire", "neut",   
                      "lymp", "mono","eosi", "baso")

#####set of covariates considered ####
covset <- cbind(covset, lneut = log(covset$neut), llymp = log(covset$lymp), 
                lmono = log(covset$mono), leosi = log(covset$eosi), 
                lbaso = log(covset$baso)) 
covset$Line <- as.factor(covset$Line)
covset$Diet <- as.factor(covset$Diet)
covset$Block <- as.factor(covset$Block)
covset$Blockorder <- as.factor(covset$Blockorder)
covset$iddam <- as.factor(covset$iddam)
covset$idsire <- as.factor(covset$idsire)
levels(covset$idsire) <- 1:11
levels(covset$iddam) <- 1:20

covset[, c("iddam", "idsire")]
detach(covset)
attach(covset)

###List of models function ####
### Case 1: no cbc data ####

colnames(covset)
full_model <- model.matrix(~Line*Diet*RFI + Block + 
                             Blockorder + Concb + RINb + 
                             Conca + RINa + lneut + llymp+ lmono + leosi + lbaso)
colnames(full_model)
rankMatrix(full_model)
list_model <- function(full_model){
  n <- dim(full_model)[2]
  variable_name <- colnames(full_model)
  variable_name <- gsub(":", "", variable_name)
  for (i in 2:8){variable_name <- gsub(i, "", variable_name)}
  test.mat <- NULL
  design.list <- vector("list", n)
  design.list[[1]] <- full_model
  for (i in 2:n) {
    design.list[[i]] <- as.matrix(full_model[,-i])
    test.mat <- rbind(test.mat, c(1,i))
  }
  
  row.names(test.mat) <-  variable_name[-1]
  
  if (any(variable_name == "dateGD1/13/01")){
    ind_dateGD <- which(variable_name == "dateGD1/13/01")
    design.list <- vector("list", n-3)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(n-4), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 2:(ind_dateGD-1)){
      design.list[[i]] <- as.matrix(full_model[,-i])
      row.names(test.mat)[i-1] <- variable_name[i]
    }
    
    design.list[[ind_dateGD]] <- full_model[,-(ind_dateGD:(ind_dateGD+3))]
    row.names(test.mat)[ind_dateGD -1] <- "dateGD"
    
    if((ind_dateGD+1)<=(n-3)){
      for(i in ((ind_dateGD+1):(n-3))){
        design.list[[i]] <- full_model[,-(i+3)]
        row.names(test.mat)[i-1] <- variable_name[i+3]    
      }
    }
  }
  
  if (any(variable_name == "dateRNA11/14/01")){
    ind_dateRNA <- which(variable_name == "dateRNA11/14/01")
    design.list <- vector("list", n-1)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(n-2), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 2:(ind_dateRNA-1)){
      design.list[[i]] <- as.matrix(full_model[,-i])
      row.names(test.mat)[i-1] <- variable_name[i]
    }
    
    design.list[[ind_dateRNA]] <- full_model[,-(ind_dateRNA:(ind_dateRNA+1))]
    row.names(test.mat)[ind_dateRNA -1] <- "dateRNA"
    
    if ((ind_dateRNA+1)<=(n-1)) {for(i in ((ind_dateRNA+1):(n-1))){
      design.list[[i]] <- full_model[,-(i+1)]
      row.names(test.mat)[i-1] <- variable_name[i+1]    
    }
    }
  }
  if (n ==2) design.list[[2]] <- rep(1, nrow(full_model))
  return(list(design.list = design.list, test.mat = test.mat))
}


## Function do all the things with input Full model

fit_model <- function(full_model, model_th){
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, 
                log.offset = log.offset, print.progress=FALSE,
                Model = "NegBin")
  result<- QL.results(fit, Plot = FALSE)
  k <- nrow(test.mat)
  name_model <- NULL
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
  dir.create(model_dir, showWarnings = FALSE)
  save(result, file = paste(model_dir,"/Model",model_th, "_result.RData", sep =""))
  save(fit, file = paste(model_dir,"/Model",model_th, "_fit.RData", sep =""))
  for(i in 1:(nrow(test.mat))){
    postscript(paste(model_dir,"/Model", 
                     model_th, row.names(test.mat)[i],".eps", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
    
    pdf(paste(model_dir,"/Model", 
              model_th, row.names(test.mat)[i],".pdf", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
  }
  print(paste("Model", model_th, sep = " "))
}


m <- 1
full_model <- model.matrix(~block+ blockorder)
full_model <- model.matrix(~idsire + iddam)
rankMatrix(full_model)

#colnames(full_model)
# list_model(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1
