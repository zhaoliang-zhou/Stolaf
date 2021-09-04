# Mercury Example

# Bayesian approach
halflife <- c(52,69,73,88,87,56)
sigma <- 8
priormu <- 80
priornu <- 8
xbar <- mean(halflife)
n <- length(halflife)
wt <- sigma^2 / (sigma^2 + n*priornu^2)
postmu <- wt*priormu + (1-wt)*xbar
postnu <- sqrt( (sigma^2)*(priornu^2) / (sigma^2 + n*priornu^2) )
theta <- seq(60,90,length=1000)
dtheta <- theta[2] - theta[1]  # approx 30/1000
yprior <- dnorm(theta,priormu,priornu)
ypost <- dnorm(theta,postmu,postnu)

# Check numerical stuff
sum(ypost*dtheta)
postnu

# Posterior mean = E(theta | data) where our model says X~N(theta,sigma)
postmu                   # exact (derived) posterior mean
sum(theta*ypost*dtheta)  # approximate (numerical) posterior mean

# Pr(theta < 80 | data)
sum(ypost[theta<80]*dtheta)  # exact (derived) probability
pnorm(80,postmu,postnu)      # approximate (numerical) probability

# 95% credible interval for theta
qnorm(.025,postmu,postnu)    # exact (derived) interval
qnorm(.975,postmu,postnu)
lowertail <- max(theta[cumsum(ypost*dtheta)<.025])  # approx (numerical)
uppertail <- min(theta[cumsum(ypost*dtheta)>.975])  #   interval
lowertail; uppertail

sumdata <- sum(halflife)
sumdatasq <- sum(halflife^2)
ylikeli <- exp((-sumdatasq + 2*theta*sumdata - n*theta^2) / (2*sigma^2))
ylikelihood <- ylikeli/(sum(ylikeli)*dtheta)  # Divide by area under curve
bayes2 <- data.frame(theta,yprior,ylikelihood,ypost)
plot(theta,ypost,type="l", ylim=c(0,max(ypost,ylikelihood,yprior)),
     xlab="Normal mean theta",
     ylab="Probability")
title("Mercury Pollution")
legend(78,.13,c("Posterior","Prior","Likelihood"),lty=c(1,2,3),cex=0.5)
lines(theta,ylikelihood,lty=3)
lines(theta,yprior,lty=2)

# Frequentist (classical) approach
z <- (xbar-80) / (sigma / sqrt(n))
pvalue <- pnorm(z,0,1)
lower.bound <- xbar - 1.96*sigma/sqrt(n)
upper.bound <- xbar + 1.96*sigma/sqrt(n)
z; pvalue; lower.bound; upper.bound
lowertail; uppertail
