# Bootstrapping - Section 11.5

# Read in data from Chance and Rossman
temp1<-read.csv("~/Stats 322 S20/Class/Data/gettysburgsample.csv")
temp2<-read.csv("~/Stats 322 S20/Class/Data/gettysburgpop.csv")
temp3<-read.csv("~/Stats 322 S20/Class/Data/heroin.csv")

library(boot)

# Investigation 4.5.1 - Gettysburg Address

# Population parameters
letters <- temp2[,2]
mu <- mean(letters)
truesd <- sd(letters)/sqrt(10)

# Sampling distribution of xbar (applies when know population)
#   Match plots and summary stats on p.366 of Chance Rossman
xbars <- rep(0,1000)
sds <- rep(0,1000)
for (i in 1:1000)  {
  samp10 <- sample(letters,10)
  xbars[i] <- mean(samp10)
  sds[i] <- sd(samp10)  }
hist(xbars)
mean(xbars)
sd(xbars)
qqnorm(xbars)
# Sampling distribution of t-statistics
tstats <- (xbars - mu) / (sds/sqrt(10))
hist(tstats)
qqplot(tstats, rt(1000,df=9))

# Original sample of size 10 from population  (Part [c], p.367)
origsamp <- temp1[,2]
thetahat <- mean(origsamp)
sdhat <- sd(origsamp)
t.test(origsamp)

# Bootstrap samples from original sample
bootsamps <- matrix(0,nrow=1000,ncol=10)
for (i in 1:1000)  {
  bootsamps[i,] <- sample(origsamp,size=10,replace=T)  }
bootmeans <- apply(bootsamps,1,mean)

# Investigate individual bootstrap samples  (Parts [d]-[e], p.367)
par(mfrow=c(2,1))
hist(bootsamps[1,])
hist(bootsamps[2,])
mean(bootsamps[1,])
mean(bootsamps[2,])

# Investigate all bootstrap samples (Parts [f]-[g], p.368)
par(mfrow=c(1,1))
hist(bootmeans)
mean(bootmeans)
sd(bootmeans)

# Rough 95% CI for mu  (Part [h], p.369)
#   Okay if bootstrap distribution looks normal
mean(bootmeans) - qt(.975,9)*sd(bootmeans)
mean(bootmeans) + qt(.975,9)*sd(bootmeans)

# Percentile bootstrap 95% CI for mu  (Described on p.369-371)
thetastar025 <- sort(bootmeans)[.025*1000]
thetastar975 <- sort(bootmeans)[.975*1000]
lb <- 2*thetahat - thetastar975
ub <- 2*thetahat - thetastar025
lb
ub

# Bootstrap inference using boot() function in R
meanboot <- function(x,i) {mean(x[i])}
boot1.out <- boot(origsamp,meanboot,R=1000)
boot.ci(boot1.out)


# Investigation 4.5.3 - Heroin Addiction

# Describe data in original sample
time <- temp3[,4]
hist(time)
thetahat <- median(time)
thetahat
IQR(time)
n <- length(time)

# Bootstrap samples from original sample
bootsamps <- matrix(0,nrow=1000,ncol=n)
for (i in 1:1000)  {
  bootsamps[i,] <- sample(time,size=n,replace=T)  }
bootmeds <- apply(bootsamps,1,median)

# Investigate individual bootstrap samples
par(mfrow=c(2,1))
hist(bootsamps[1,])
hist(bootsamps[2,])
mean(bootsamps[1,])
mean(bootsamps[2,])
sd(bootsamps[1,])
sd(bootsamps[2,])

# Investigate all bootstrap samples
par(mfrow=c(1,1))
hist(bootmeds)
table(bootmeds)
mean(bootmeds)
sd(bootmeds)

# Rough 95% CI for true median
mean(bootmeds) - qt(.975,n-1)*sd(bootmeds)
mean(bootmeds) + qt(.975,n-1)*sd(bootmeds)

# Percentile bootstrap 95% CI for true median
thetastar025 <- sort(bootmeds)[.025*1000]
thetastar975 <- sort(bootmeds)[.975*1000]
lb <- 2*thetahat - thetastar975
ub <- 2*thetahat - thetastar025
lb
ub

# Bootstrap inference using boot() function in R
medboot <- function(x,i) {median(x[i])}
boot1.out <- boot(time,medboot,R=1000)
boot.ci(boot1.out)


# Inference about trimmed mean
thetahat <- mean(time,trim=.25)
thetahat
boottrmeans <- apply(bootsamps,1,mean,trim=.25)

# Investigate all bootstrap samples
par(mfrow=c(1,1))
hist(boottrmeans)
mean(boottrmeans)
sd(boottrmeans)

# Rough 95% CI for true trimmed mean
mean(boottrmeans) - qt(.975,n-1)*sd(boottrmeans)
mean(boottrmeans) + qt(.975,n-1)*sd(boottrmeans)

# Percentile bootstrap 95% CI for true trimmed mean
thetastar025 <- sort(boottrmeans)[.025*1000]
thetastar975 <- sort(boottrmeans)[.975*1000]
lb <- 2*thetahat - thetastar975
ub <- 2*thetahat - thetastar025
lb
ub

# Bootstrap inference using boot() function in R
trmeanboot <- function(x,i) {mean(x[i],trim=.25)}
boot1.out <- boot(time,trmeanboot,R=1000)
boot.ci(boot1.out)
