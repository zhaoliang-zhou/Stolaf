# Nursing Homes in New Mexico - Example from Section 8.6
# MCMC Gibbs Sampling - Section 12.5 (Examples 12.5.2, 12.5.3)

# The data - medical in-patient days (in hundreds)
x <- c(128,281,291,238,155,148,154,232,316,96,146,151,100,213,208,157,48,217)
n <- length(x)
xbar <- mean(x)
snsq <- (n-1)*var(x)

# Prior hyperparameters (read derivation process on p.500)
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
m <- 100  # burn-in period
m0 <- 1000  # iterations per chain (after burn-in)
mustart <- c(182.17, 227, 272, 137, 82)  # k=5 initial values
# Note: choice of starting values described in Example 12.5.2

muall <- matrix(0,ncol=5,nrow=m0)
tauall <- matrix(0,ncol=5,nrow=m0)
muburn <- rep(0,length=5*m)
tauburn <- rep(0,length=5*m)
for(k in 1:5)  {
  mu <- rep(0,length=m+m0)
  tau <- rep(0,length=m+m0)
  mu[1] <- mustart[k]
  tau[1] <- rgamma(1,alpha1+0.5,beta1+0.5*(lambda1*(mu[1]-mu1)^2))
  for(i in 2:(m+m0))  {
    mu[i] <- rnorm(1,mu1,sqrt(1/(tau[i-1]*lambda1)))
    tau[i] <- rgamma(1,alpha1+0.5,beta1+0.5*(lambda1*(mu[i]-mu1)^2))  }
  muburn[((k-1)*m+1):(k*m)] <- mu[1:m]  # save burn-in to check convergence
  tauburn[((k-1)*m+1):(k*m)] <- tau[1:m]
  muall[,k] <- mu[-(1:m)]  # remove burn-in and place in kth column
  tauall[,k] <- tau[-(1:m)]  }

# Plot progression of Markov chain 5
par(mfrow=c(2,2))
plot(1:100,muburn[401:500],type="l")  # just burn-in
plot(1:1100,c(muburn[401:500],muall[,5]),type="l")
# Plot progression of Markov chain 1
plot(1:100,muburn[1:100],type="l")  # just burn-in
plot(1:1100,c(muburn[1:100],muall[,1]),type="l")
par(mfrow=c(1,1))



# F statistic to determine if burn-in can be stopped
chain <- as.factor(rep(1:5,each=m))
anova(lm(muburn~chain))

anova(lm(tauburn~chain))

# Simulated values of mu and tau
mu <- as.vector(muall)
tau <- as.vector(tauall)
U <- sqrt((lambda1*alpha1)/beta1)*(mu-mu1)

# Emulate quantile plots in Figure 12.7
par(mfrow=c(1,2))
qqplot(U,rt(5*m0,2*alpha1))
qqplot(tau,rgamma(5*m0,alpha1,beta1))
par(mfrow=c(1,1))

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


mean(mu+1.645/sqrt(tau))  # Z = 0.95 quantile of original normal distribution

mu1+1.645*sqrt(beta1)*gamma(alpha1-0.5)/gamma(alpha1)  # actual


# Compute simulation standard error of Z
Z5 <- rep(0,length=5)
for(k in 1:5)  {
  Z5[k] <- mean(muall[,k]+1.645/sqrt(tauall[,k]))  }
S <- sqrt(4*var(Z5)/5)
sigmahat <- sqrt(m0*S^2)
simSE <- sigmahat/sqrt(5*m0)

# Compare crude 95% prediction interval from MCMC to actual
lb <- mu1 - sqrt(beta1/(lambda1*alpha1))*qt(.975,2*alpha1)
ub <- mu1 + sqrt(beta1/(lambda1*alpha1))*qt(.975,2*alpha1)

sort(mu)[floor(.025*5*m0)]
sort(mu)[ceiling(.975*5*m0)]


