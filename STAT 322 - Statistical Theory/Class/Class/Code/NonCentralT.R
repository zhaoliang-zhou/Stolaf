# Check out the noncentral t distribution

psi <- 6
n <- 10
num <- rnorm(10000,psi,1)
denom <- sqrt(rchisq(1000,n-1)/(n-1))
u <- num/denom

nont <- rt(10000,df=n-1,ncp=psi)
tpluspsi <- rt(10000,df=n-1) + psi
summary(nont)
summary(tpluspsi)
summary(u)
sd(u)
sd(nont)
sd(tpluspsi)
n/(n-2)
plot(density(tpluspsi))
lines(density(nont),lty=2)
lines(density(u),lty=4)

