#source("Utils.R")                  # Loads utility functions (assumes it is in the same directory)

# Number of people in the populations
N <- 100000 

# Vaccine coverage
a.phi <- betaPar2(.434,.6,.95)$res1
b.phi <- betaPar2(.434,.6,.95)$res2

# Baseline probabilities of clinical outcomes
# 1. Influenza infection
# 2. GP visits
# 3. Minor complications (repeat visit)
# 4. Major complications (pneumonia)
# 5. Hospitalisations
# 6. Death
# 7. Adverse events due to vaccination
N.outcomes <- 7
mu.beta <- c(.0655,.273,.401,.0128,.00038,.00075)
upp.beta <- c(.111,.51,.415,.0197,.00067,.000132)
sd.beta <- c(NA,NA,NA,NA,.0001,.00028)
a.beta <- b.beta <- numeric()
for (i in 1:4) {
    a.beta[i] <- betaPar2(mu.beta[i],upp.beta[i],.975)$res1
    b.beta[i] <- betaPar2(mu.beta[i],upp.beta[i],.975)$res2
}
for (i in 5:6) {
    a.beta[i] <- lognPar(mu.beta[i],sd.beta[i])$mulog
    b.beta[i] <- 1/lognPar(mu.beta[i],sd.beta[i])$sigmalog^2
}
a.beta[N.outcomes] <- betaPar(.1,.05)$a
b.beta[N.outcomes] <- betaPar(.1,.05)$b

# Decrease in risk of infection due to vaccination
mu.rho <- lognPar(.69,.05)$mulog; tau.rho <- 1/lognPar(.69,.05)$sigmalog^2

# Treatment with antibiotics after GP visit
mu.gamma <- c(.42,.814);  sd.gamma <- c(.0015,.004)
a.gamma <- b.gamma <- numeric()
for (i in 1:2) {
    a.gamma[i] <- betaPar(mu.gamma[i],sd.gamma[i])$a
    b.gamma[i] <- betaPar(mu.gamma[i],sd.gamma[i])$b
}

# Number of prescriptions of antibiotics
a.delta <- 7

# Taking OTC
a.xi <- betaPar(.95,.005)$a
b.xi <- betaPar(.95,.005)$b

# Being off work
a.eta <- betaPar(.9,.005)$a
b.eta <- betaPar(.9,.005)$b

# Length of absence from work for influenza
mu.lambda <- lognPar(2.9,1.25)$mulog
tau.lambda <- 1/lognPar(2.9,1.25)$sigmalog^2

# Costs (N.resources = 8)
# 1. Cost of GP visit
# 2. Cost of hospital episode
# 3. Cost of vaccination
# 4. Cost of time off for individuals to receive vaccination
# 5. Cost of days work absence due to influenza
# 6. Cost of antibiotics
# 7. Cost of OTC treatments
# 8. Cost of travel to receive vaccination
N.resources <- 8
m.psi <- c(20.66,2656,7.24,10.16,46.27,3.81,1.6,.81)
sd.psi <- c(5.015,440.75,1.81,2.54,11.57,.955,.4,.2)
sd.psi <- .25*m.psi
mu.psi <- tau.psi <- rep(0,N.resources)
for (i in 1:N.resources) {
    mu.psi[i] <- lognPar(m.psi[i],sd.psi[i])$mulog
    tau.psi[i] <- 1/lognPar(m.psi[i],sd.psi[i])$sigmalog^2
}

# Quality of life weights (N.outcomes = 7)
# 1. Influenza infection
# 2. GP visits (no QoL loss)
# 3. Minor complications (repeat visit, no QoL loss)
# 4. Major complications (pneumonia)
# 5. Hospitalisations
# 6. Death
# 7. Adverse events due to vaccination
m.omega <- c(4.27,0,0,6.35,6.35,15.29,.55)
sd.omega <- c(1.38,0,0,1.5875,1.5875,3.6,.15)
mu.omega <- tau.omega <- rep(0,N.outcomes)
for (i in c(1,4,5,6,7)) {
    mu.omega[i] <- lognPar(m.omega[i],sd.omega[i])$mulog
    tau.omega[i] <- 1/lognPar(m.omega[i],sd.omega[i])$sigmalog^2
}

