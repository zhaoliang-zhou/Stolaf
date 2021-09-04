# Example 2.3.7

j <- 1:11
B <- (j-1)/10   # proportion of times imipramine patient successful = Pr(E1 | Bj)
prior <- rep(1/11,11)  # prior probabilities = Pr(Bj)
priormean <- sum(B*prior) # prior mean = Pr(E1)
priormean
# this is also the predictive probability that the next patient is successful

# Assume we learn 22 of first 40 imipramine patients non-relapsers
s <- 22  # s is a success (non-relapser)
f <- 18  # f is a failure (relapser)
likelihood <- choose(s+f,s)*(B^s)*((1-B)^f)  # Pr(A | Bj)
normconst <- sum(prior*likelihood)  # Pr(A) = normalizing constant
posterior <- prior*likelihood/normconst  # Pr(Bj | A)
posterior
sum(posterior)

# Find mean of posterior distribution = Pr(E41 | A)
postmean <- sum(B*posterior)
postmean
# Find posterior probability that prob of success is above 50% (p > .5)
sum(posterior[7:11])

bayes10 <- data.frame(B,prior,likelihood,posterior)
bayes10

# Plot - superimpose prior (dashed) and posterior (solid) bar charts.
plot(rep(B[1],2),c(0,posterior[1]),xlim=c(0,1),
     ylim=c(0,max(posterior,prior)),type="l",xlab="Proportion of non-relapsers",
     ylab="Probability (solid=posterior, dashed=prior)")
title("Example 2.3.7")
for(i in 2:11)  {
  lines(rep(B[i],2),c(0,posterior[i]))
}
for(i in 1:11)  {
  lines(rep(B[i]+.02,2),c(0,prior[i]),lty=2)
}

# Example 2.3.8 - repeat 2.3.7 with new prior opinion about p
prior <- c(0,0.39,0.20,0.13,0.09,0.07,0.05,0.03,0.02,0.02,0)
priormean <- sum(B*prior)
priormean
likelihood <- choose(s+f,s)*(B^s)*((1-B)^f)  # Pr(A | Bj)
normconst <- sum(prior*likelihood)  # Pr(A)
posterior <- prior*likelihood/sum(prior*likelihood)  # Pr(Bj | A)
posterior

# Find mean of posterior distribution = Pr(E41 | A)
postmean <- sum(B*posterior)
postmean
# Find posterior probability that prob of success is above 50% (p > .5)
sum(posterior[7:11])

bayes11 <- data.frame(B,prior,likelihood,posterior)
bayes11

# Plot - superimpose prior (dashed) and posterior (solid) bar charts.
plot(rep(B[1],2),c(0,posterior[1]),xlim=c(0,1),
     ylim=c(0,max(posterior,prior)),type="l",xlab="Proportion of non-relapsers",
     ylab="Probability (solid=posterior, dashed=prior)")
title("Example 2.3.8")
for(i in 2:11)  {
  lines(rep(B[i],2),c(0,posterior[i]))
}
for(i in 1:11)  {
  lines(rep(B[i]+.02,2),c(0,prior[i]),lty=2)
}

# Plot two priors from Example 2.3.7 and 2.3.8
plot(rep(B[1],2),c(0,bayes10[1,2]),xlim=c(0,1),
     ylim=c(0,max(bayes10[,2],bayes11[,2])),type="l",
     xlab="Prior proportion of non-relapsers",
     ylab="Probability (solid=Ex2.3.7, dashed=Ex2.3.8)")
title("Priors from two examples")
for(i in 2:11)  {
  lines(rep(B[i],2),c(0,bayes10[i,2]))
}
for(i in 1:11)  {
  lines(rep(B[i]+.02,2),c(0,bayes11[i,2]),lty=2)
}

# Plot two posteriors from Example 2.3.7 and 2.3.8
plot(rep(B[1],2),c(0,bayes10[1,4]),xlim=c(0,1),
     ylim=c(0,max(bayes10[,4],bayes11[,4])),type="l",
     xlab="Posterior proportion of non-relapsers",
     ylab="Probability (solid=Ex2.3.7, dashed=Ex2.3.8)")
title("Posteriors from two examples")
for(i in 2:11)  {
  lines(rep(B[i],2),c(0,bayes10[i,4]))
}
for(i in 1:11)  {
  lines(rep(B[i]+.02,2),c(0,bayes11[i,4]),lty=2)
}