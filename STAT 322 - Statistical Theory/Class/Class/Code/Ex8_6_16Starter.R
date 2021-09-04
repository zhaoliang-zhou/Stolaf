# Problem 8.6.16

# The data - lactic acid concentrations in cheese
x10 <- c(0.86,1.53,1.57,1.81,0.99,1.09,1.29,1.78,1.29,1.58)
x20 <- c(1.68,1.9,1.06,1.3,1.52,1.74,1.16,1.49,1.63,1.99,
         1.15,1.33,1.44,2.01,1.31,1.46,1.72,1.25,1.08,1.25)
x <- c(x10,x20)
n <- length(x); n
xbar <- mean(x); xbar
snsq <- (n-1)*var(x); snsq

# Prior hyperparameters from Example 8.6.2
mu0 <- 1
lambda0 <- 1
alpha0 <- 0.5
beta0 <- 0.5

# Marginal prior for tau
x1 <- seq(0.1,20,length=1000)
y1 <- dgamma(x1,alpha0,beta0)
plot(x1,y1,type="l",xlab="tau",ylab="density",
  main="Marginal prior distribution of tau")

# Marginal prior for mu
x2 <- seq(-3,5,length=1000)
dx2 <- x2[2]-x2[1]
U <- sqrt((lambda0*alpha0)/beta0)*(x2-mu0)
y2 <- dt(U,2*alpha0)
plot(x2,y2/(sum(y2)*dx2),type="l",xlab="mu",ylab="density",
  main="Marginal prior distribution of mu")

# Joint prior distribution of mu and tau
x3 <- seq(.1,20,length=100)
y3 <- seq(-3,5,length=100)
xy <- expand.grid(x3,y3)
z3 <- dgamma(xy[,1],alpha0,beta0)*dnorm(xy[,2],mu0,1/sqrt(lambda0*xy[,1]))
z3mat <- matrix(z3,ncol=100)
persp(x3,y3,z3mat,xlab="tau",ylab="mu",zlab="joint density")
contour(x3,y3,z3mat,xlab="tau",ylab="mu")
image(x3,y3,z3mat,xlab="tau",ylab="mu")

# Marginal prior for mu based on joint distribution above
marg.mu <- apply(z3mat,2,sum)
#sum across columns of grid/matrix
dy3 <- y3[2]-y3[1]
ht <- marg.mu/(sum(marg.mu)*dy3)
plot(x2,y2/(sum(y2)*dx2),type="l",xlab="mu",ylab="density",
  main="Marginal prior distribution of mu (collapsed dotted)")
lines(y3,ht,lty=2)

# Prior probability interval for mu
lb <- mu0 - sqrt(beta0/(lambda0*alpha0))*qt(.975,2*alpha0)
ub <- mu0 + sqrt(beta0/(lambda0*alpha0))*qt(.975,2*alpha0)
lb
ub


# Posterior hyperparameters
mu1 <- (lambda0*mu0+n*xbar)/(lambda0+n); mu1
lambda1 <- lambda0 + n ; lambda1
alpha1 <- alpha0 + (n/2) ; alpha1
beta1 <- beta0 + (snsq/2) + ((n*lambda0*(xbar-mu0)^2)/(2*(lambda0+n))); beta1

# Marginal posterior for tau
x4 <- seq(0.1,20,length=1000)
y4 <- dgamma(x4,alpha1,beta1)
plot(x4,y4,type="l",xlab="tau",ylab="density",
  main="Marginal posterior distribution of tau (prior dotted)")
lines(x1,y1,lty=2)
 
# Marginal posterior for mu
x5 <- seq(-1,3,length=1000)
dx5 <- x5[2]-x5[1]
U <- sqrt((lambda1*alpha1)/beta1)*(x5-mu1)
y5 <- dt(U,2*alpha1)
#y5 is the marginal prior based on U
plot(x5,y5/(sum(y5)*dx5),type="l",xlab="mu",ylab="density",
  main="Marginal posterior distribution of mu (prior dotted)")
lines(x2,y2/(sum(y2)*dx2),lty=2)

# Joint posterior distribution of mu and tau
x6 <- seq(.1,20,length=100)
y6 <- seq(-3,5,length=100)
xy <- expand.grid(x6,y6)
z6 <- dgamma(xy[,1],alpha1,beta1)*dnorm(xy[,2],mu1,1/sqrt(lambda1*xy[,1]))
z6mat <- matrix(z6,ncol=100)
persp(x6,y6,z6mat,xlab="tau",ylab="mu",zlab="joint density")
contour(x6,y6,z6mat,xlab="tau",ylab="mu",
  main="Joint distributions (prior=dotted, posterior=solid)")
     priorcontour <- contourLines(x3,y3,z3mat)
     templines <- function(clines) {
       lines(clines[[2]], clines[[3]],lty=2)
     }  
     invisible(lapply(priorcontour, templines))
image(x6,y6,z6mat,xlab="tau",ylab="mu")

# Marginal posterior for mu based on joint distribution above
marg.mu <- apply(z6mat,2,sum)
dy6 <- y6[2]-y6[1]
ht <- marg.mu/(sum(marg.mu)*dy6)
plot(x5,y5/(sum(y5)*dx5),type="l",xlab="mu",ylab="density",
  main="Marginal distribution of mu (prior and 2 posteriors)")
lines(y6,ht,lty=2)
lines(x2,y2/(sum(y2)*dx2),lty=4)
#y2 is marginal prior, y5 is marginal posterior, ht is marginal based on joint)

# Posterior probability interval for mu (90% interval)
lb <- mu1 - sqrt(beta1/(lambda1*alpha1))*qt(.95,2*alpha1)
ub <- mu1 + sqrt(beta1/(lambda1*alpha1))*qt(.95,2*alpha1)
lb
ub

# Posterior probability interval for mu based on marginal posterior
lbnum <- sum(cumsum(ht)/sum(ht) < .05)
lbemp <- y6[lbnum]
ubnum <- sum(cumsum(ht)/sum(ht) < .95)
ubemp <- y6[ubnum+1]
lbnum
lbemp
ubnum
ubemp

# Compare with 95% confidence interval for mu
t.test(x,conf.level=.90)$conf.int
