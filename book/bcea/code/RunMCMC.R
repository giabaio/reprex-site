# Loads the package to run OpenBUGS or JAGS from R 
library(R2OpenBUGS)
library(R2jags)   
# Defines the current as the working directory
working.dir <- paste(getwd(),"/",sep="")
# Launches the file Utils.R which contains useful functions used throughout this script
source("https://gianluca.statistica.it/book/bcea/code/Utils.R")
# Loads the data into R (assumes the file is stored in the working directory - if not the full path can be provided)
source("https://gianluca.statistica.it/book/bcea/code/LoadData.R")                      

# Defines the data list to be passed to BUGS/JAGS
data <- list(
  "N","a.phi","b.phi","mu.rho","tau.rho","a.beta","b.beta","a.gamma","b.gamma","mu.omega",
  "tau.omega","mu.psi","tau.psi","N.outcomes","N.resources","mu.lambda","tau.lambda",
  "a.xi","b.xi","a.eta","b.eta","a.delta"
) 
    
# Defines the file with the model code
filein <- "vaccine.txt"

# Defines the quantities to be monitored (stored)
params <- c(
  "beta","phi","omega","rho","Infected","GP","Repeat.GP","Pneumonia","Hospital","Death",
  "Mild.Compl","Trt","Adverse.events","n","gamma","delta","psi","lambda","pi","xi","eta"
)

# Generates the initial values
inits <- function(){
  list(
    phi=runif(1),beta=runif(N.outcomes,0,1),rho=c(NA,runif(1)),gamma=runif(2,0,1),delta=rpois(1,2),
    omega=c(runif(1),NA,NA,runif(1),NA,runif(2,0,1)),psi=runif(N.resources,0,10),lambda=runif(1),
    eta=runif(1),xi=runif(1)
  )
}

# Defines the number of iteration, burn-in and thinning, and runs BUGS or JAGS
n.iter <- 100000
n.burnin <- 9500
n.thin <- floor((n.iter-n.burnin)/500)

# 1. This runs OpenBUGS
vaccine <- bugs(data, inits, params, model.file=filein,n.chains=2, n.iter, n.burnin, n.thin, DIC=FALSE, working.directory=working.dir)

# 2. This runs JAGS
vaccine <- jags(data, inits, params, model.file=filein,n.chains=2, n.iter, n.burnin, n.thin, DIC=FALSE, working.directory=working.dir, progress.bar="text")
    
# Prints the summary stats and attaches the results to the R workspace
print(vaccine,digits=3,intervals=c(0.025, 0.975))

# In OpenBUGS:
attach.bugs(vaccine)
# In JAGS:
attach.jags(vaccine)
