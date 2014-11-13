source("simulation_loadModel7.R")
source("simulation_criteria.R")
source("simulation_listmodel.R")
source("simulation_fitmodel.R")
source("simulation_simcount.R")
source("simulation_runnrep.R")
## simulation replication #####
pm1 <- proc.time()
size <- 5000
# vnrep <- 1:10
#  vnrep <- 11:20
# vnrep <- 21:30
# vnrep <- 31:40
# vnrep <- 41:50
# vnrep <- 51:60
# vnrep <- 61:70
# vnrep <- 71:80
# vnrep <- 81:90
vnrep <- 91:100

criteria <- 1 #pvalue05
for(nrep in vnrep) # nrep <- 500
{
  simdat <- simcount(nrep, size)
  save(simdat, file = paste0("simdat/simdat_", nrep, ".RData"))
  print(paste0("nrep = ", nrep ))
  runnrep(simdat, criteria)
}


proc.time()-pm1