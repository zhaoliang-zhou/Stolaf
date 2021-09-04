# Example 12.6.4 from DeGroot and Schervish p843
table1033 <- c(-2.46,-2.11,-1.23,-0.99,-0.42,-0.39,-0.21,-0.15,-0.10,
              -0.07,-0.02,0.27,0.40,0.42,0.44,0.70,0.81,0.88,1.07,1.39,1.40,
              1.47,1.62,1.64,1.76)
hist(table1033)
thetahat <- IQR(table1033);thetahat
n <- length(table1033);n

#using IQR function is a little different than the book
sort_tab<-sort(table1033)
bookIQR=sort_tab[19]-sort_tab[6];bookIQR

# Bootstrap samples from original sample
bootsamps <- matrix(0,nrow=10000,ncol=n)
for (i in 1:10000)  {
  bootsamps[i,] <- sample(table1033,size=n,replace=T)  }
bootiqr <- apply(bootsamps,1,IQR)

# Investigate all bootstrap samples
par(mfrow=c(1,1))
hist(bootiqr)
mean(bootiqr)
sd(bootiqr)

# Percentile bootstrap 90% CI for IQR
a <- sort(bootiqr/thetahat)[.05*10000]
b <- sort(bootiqr/thetahat)[.95*10000]
lb <- thetahat/b
ub <- thetahat/a
lb
ub

# Bootstrap inference using boot() function in R
iqrboot <- function(x,i) {IQR(x[i])}
boot1.out <- boot(table1033,iqrboot,R=10000)
boot.ci(boot1.out, conf=.90)


# Example 12.6.6 from DeGroot and Schervish
lactic <- c(0.86,1.53,1.57,1.81,0.99,1.09,1.29,1.78,1.29,1.58)
hist(lactic)
m <- median(lactic);m
y <- mad(lactic,constant=1);y
n <- length(lactic);n

# Bootstrap samples from original sample
bootsamps <- matrix(0,nrow=10000,ncol=n)
for (i in 1:10000)  {
  bootsamps[i,] <- sample(lactic,size=n,replace=T)  }
mstar <- apply(bootsamps,1,median)
ystar <- apply(bootsamps,1,mad,constant=1)

# Investigate all bootstrap samples
par(mfrow=c(1,1))
hist(mstar)
mean(mstar)
sd(mstar)
par(mfrow=c(1,1))
hist(ystar)
mean(ystar)
sd(ystar)

# Percentile-t bootstrap 90% CI for median
d2 <- -sort((mstar-m)/ystar)[.05*10000]
d1 <- sort((mstar-m)/ystar)[.95*10000]
lb <- m-d1*y
ub <- m+d2*y
lb
ub

# Percentile bootstrap 90% CI for median
d2 <- -sort(mstar-m)[.05*10000]
d1 <- sort(mstar-m)[.95*10000]
lb <- m-d1
ub <- m+d2
lb
ub

# Bootstrap inference using boot() function in R
medboot <- function(x,i) {median(x[i])}
boot1.out <- boot(lactic,medboot,R=10000)
boot.ci(boot1.out,conf=.90)
