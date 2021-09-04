# Nursing Homes in New Mexico - Example from Section 8.6
# MCMC Gibbs Sampling - Section 12.5 (Examples 12.5.2, 12.5.3)

# The data - medical in-patient days (in hundreds)
x <- c(128,281,291,238,155,148,154,232,316,96,146,151,100,213,208,157,48,217)
n <- length(x)
xbar <- mean(x)
snsq <- (n-1)*var(x)

# Prior hyperparameters
mu0 <- 200
lambda0 <- 2
alpha0 <- 2
beta0 <- 6300

# Posterior hyperparameters
mu1 <- (lambda0*mu0+n*xbar)/(lambda0+n)
lambda1 <- lambda0 + n
alpha1 <- alpha0 + (n/2)
beta1 <- beta0 + (snsq/2) + ((n*lambda0*(xbar-mu0)^2)/(2*(lambda0+n)))

# Gibbs sampling algorithm
m <- 1000  # chain length
mu <- rep(0,length=m)
tau <- rep(0,length=m)

##
## Your turn
## Fill in vectors of mus and taus according to Gibbs sampling algorithm!
##

#Initial values
mu[1] <- #STARTINGVALUE
tau[1] <- #SIMULATED FROM CONDITIONAL
  
for(i in 2:m)  {
   }


# Plot progression of Markov chains
par(mfrow=c(2,1))
plot(1:m,mu,type="l")
plot(1:m,tau,type="l")

# Sketch posterior marginal distributions of mu and tau
x4 <- seq(0.000001,.001,length=1000)
y4 <- dgamma(x4,alpha1,beta1)
plot(x4,y4,type="l",xlab="tau",ylab="density",
  main="Marginal posterior distribution of tau (MCMC dotted)")
lines(density(tau),lty=2)

x5 <- seq(0,400,length=1000)
dx5 <- x5[2]-x5[1]
U <- sqrt((lambda1*alpha1)/beta1)*(x5-mu1)
y5 <- dt(U,2*alpha1)
plot(x5,y5/(sum(y5)*dx5),type="l",xlab="mu",ylab="density",
  main="Marginal posterior distribution of mu (MCMC dotted)")
lines(density(mu),lty=2)

# Use MCMC results to estimate quantities of interest
mean(mu)  
mu1
mean(tau)
alpha1/beta1

# Compare crude 95% prediction interval from MCMC to actual
lb <- mu1 - sqrt(beta1/(lambda1*alpha1))*qt(.975,2*alpha1)
lb
ub <- mu1 + sqrt(beta1/(lambda1*alpha1))*qt(.975,2*alpha1)
ub
sort(mu)[floor(.025*m)]
sort(mu)[ceiling(.975*m)]
