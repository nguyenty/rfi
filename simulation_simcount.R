simcount <- function(nrep, size){
  set.seed(nrep)
  J <- dim(used_beta)[1]
  #used_gene is the list of indexes of genes from the original data 12280 
  s <- sample( dim(used_beta)[1], size = size)
  s <- s[order(s)]
  s_mu <- mu[s,]
  s_omega <- used_omega[s]
  s_degene <- intersect(degene, s) 
  #length(s_degene)/length(s)
  
  ##Sim counts data
  
  y <- array(0, dim = c(size,31))
  
  # nrep <- 100;size <- 5000 ; j <-5000; k <- 1
  
  for(j in 1:size){ # j <- 1; k <- 1
    repeat{
      for(k in 1:31){
        #set.seed((32*nrep+k)*size + j)
        y[j,k]  <- rnbinom(n=1, size=1/s_omega[j], mu=s_mu[j,k])
      }
      if (mean(y[j,])>8& sum(y[j,]>0)>3) break
    }
  }
  
  
  
  # nrep <- 1
  
  return(list(s = s,  s_degene = s_degene, 
                        y = y, nrep = nrep))
}
