# Example 5.8.5 = Example 2.3.7 with continuous distributions
s=22
f=18
# Doctors prior
priora <- 1.18
priorb <- 2.35
posta <- priora + s
postb <- priorb + f
x <- seq(0,1,length=1000)
yprior <- dbeta(x,priora,priorb)
ypost <- dbeta(x,posta,postb)
ylikeli <- dbinom(s,s+f,x)
ylikelihood <- ylikeli/(sum(ylikeli*.001))  # Divide by area under curve (approx)
plot(x,ypost,type="l",xlab="Proportion of non-relapsers",
     ylab="Probability (solid=posterior, dash=prior, dots=likelihood)")
lines(x,ylikelihood,lty=3)
lines(x,yprior,lty=2)

# Uniform prior
priora <- 1
priorb <- 1
posta <- priora + s
postb <- priorb + f
x <- seq(0,1,length=1000)
yprior <- dbeta(x,priora,priorb)
ypost <- dbeta(x,posta,postb)
ylikeli <- dbinom(s,s+f,x)
ylikelihood <- ylikeli/(sum(ylikeli*.001))
plot(x,ypost,type="l",xlab="Proportion of non-relapsers",
     ylab="Probability (solid=posterior, dash=prior, dots=likelihood)")
lines(x,ylikelihood,lty=3)
lines(x,yprior,lty=2)
