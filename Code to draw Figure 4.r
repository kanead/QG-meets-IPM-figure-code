# In this bit of code I run a dynamic version of the Breeders equation
# In (a) and (d) there is no density-dependence in E. The population grows
# hyper-exponentially (d) and the phenotype and its components evolve 
# relativel slowly
# in (b) and (e) I incldue density-dependent selection. This increases 
# the rate of evolution but slows the population growth rate and increase in
# population size.
# in (c) and (f) adaptive plasticity is included in the development function for E
# this slows the rate of evolution as the phenotype can converge on its optimal value
# faster via plasticity than when plasticity is not operating.

# in the remaining three models we examine how populations respond to 
# a perturbation when no plasticity, adaptive plasticity and non-adaptive plasticity is operating.

# a little bit of book keeping
graphics.off()
rm(list=ls())
# a function to calculate the first four central, and the first
# four moments of a distribution
moments.fun <- function(x,z){
  mu <- sum(z*x)/sum(x)
  va <- sum(z*z*x)/sum(x)-mu^2
  sk <- (sum(z*z*z*x)/sum(x)-3*mu*va-mu^3)/(sqrt(va)^3)
  m4 <- sum(z*z*z*z*x)/sum(x)
  m3 <- sum(z*z*z*x)/sum(x)
  m2 <- sum(z*z*x)/sum(x)
  ku <- (m4-4*mu*m3+6*(mu^2)*m2-3*(mu^4))/(va^2)-3
  res <- c(mu,va,sk,ku,m2,m3,m4)
  names(res) <- c('Mean','Variance','Skew','Kurtosis','Moment 2','Moment 3','Moment 4')
  return(res)
}
# Define fitness function
fitness.fun <- function(x1,int,slp.z1,slp.N,popN) int+slp.z1*x1+slp.N*popN

# define Gaussian inheritance function
G.fun <- function(z, zz, mu.a, mu.b, mu.N, sigma.a, Num) {
  mu.z <- mu.a+mu.b*z+mu.N*Num
  sigma.z2 <- sigma.a
  sigma.z <- sqrt(sigma.z2)
  temp1 <- sqrt(2*pi)*sigma.z
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  return(exp(-temp2)/temp1)
}

# load required libraries (although multicore might not be use
library(mvtnorm)
library(Matrix)
# set up global parameters and choose initial values
# number of classes in the matrix
n <- 400
# values of breeding value distribution
g <- seq(1,25,length.out=n)
# create the combinations of x and y
G <- rep(g,each=n)
e <- seq(1,25,length.out=n)
E <- rep(e,times=n)
# I need to create a vector the same length of G and E with phenotypic values in it
Z <- G+E
Z.for.plot <- matrix(NA,n,n)
for (i in 1:n){
  Z.for.plot[,i] <- i:(i+(n-1))
}
range.Z <- range(Z)
G.count <- rep(1:n,each=n)
E.count <- rep(1:n,times=n)
Z.cols <- cbind(G,E)
# number of generations for the simulation
n.gens <- 100
# start with the initial distributions
mu <- c(7,12) # G,E
sigma <- matrix(c(1,0,0,1),2,2) # G,E and their covariances
Z.res.fit <- Z.res <- array(NA,c(n*n,n.gens))
# distribution of breeding values, environmental component and phenotype at time 1
Z.res[,1] <- dmvnorm(Z.cols,mean=mu,sigma=sigma)
Z.res[,1] <- Z.res[,1]/sum(Z.res[,1])

# parameters for fitness function
int <- 0.3
slp.z1 <- 0.1
slp.N <- 0
# Run the first simulation
for (i in 2:n.gens){
  popN <- sum(Z.res[,i-1]) # calculate pop size
  w <- fitness.fun(G+E,int,slp.z1,slp.N,popN) # estimate fitness
  Z.res.fit[,i-1] <- w*Z.res[,i-1] # fitness in next time step
  sizes <- tapply(Z.res.fit[,i-1],G.count,sum) # fitness of G
  Es <- tapply(Z.res.fit[,i-1],E.count,sum) # Fitness of each E
  sizes <- rep(sizes,each=n)
  dist <- dnorm(e,mu[2],sigma[4])
  dist <- dist/sum(dist)
  dist <- rep(dist,times=n)
  Z.res[,i] <- dist*sizes
}
R0 <- apply(Z.res,2,sum)
R0 <- R0[-1]/R0[-n.gens]

# start the graphics window
quartz()
par(mfrow=c(4,3),mar=c(1,4,1,1),oma=c(3,0,0.1,0.1))
# calculate moments of z, e and a
moment.res <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
for (i in 1:n.gens){
  moment.res[,i] <- moments.fun(Z.res[,i],G)[1:4]
  moment.res.z[,i] <- moments.fun(Z.res[,i],Z)[1:4]
  moment.res.e[,i] <- moments.fun(Z.res[,i],E)[1:4]
}
run.a <- Z.res
# plot the moments for the first model
plot(1:n.gens,moment.res[1,1:n.gens],type='l',xlab='',ylab='',ylim=c(0,30),axes=FALSE) #ylim=c(ranger[1],ranger[2])
axis(side=2,labels=TRUE)
axis(side=1,labels=TRUE)
box()
mtext('Mean',side=2,line=2.5)
lines(1:n.gens,moment.res.z[1,1:n.gens],col='red')
lines(1:n.gens,moment.res.e[1,1:n.gens],col='blue')
text(8,28,'(a)',cex=1.5)
par(mfg=c(2,1))
plot(1:(n.gens-1), R0[1:(n.gens-1)],type='l',xlab='',ylab='',axes=FALSE)#,ylim=c(150,250))
box()
axis(side=2,labels=TRUE)
mtext(side=2,line=2.5,'R0')
axis(side=1,labels=TRUE)
text(8,2.55,'(d)',cex=1.5)
##################################################
# now add in some dd in the selection function -- model (b). Define parameters
int <- 0.3
slp.z1 <- 0.1
slp.N <- -0.01
# run the simulation
for (i in 2:n.gens){
  popN <- sum(Z.res[,i-1]) # calculate pop size
  w <- fitness.fun(G+E,int,slp.z1,slp.N,popN) # estimate fitness
  Z.res.fit[,i-1] <- w*Z.res[,i-1] # fitness in next time step
  sizes <- tapply(Z.res.fit[,i-1],G.count,sum) # fitness of G
  Es <- tapply(Z.res.fit[,i-1],E.count,sum) # Fitness of each E
  sizes <- rep(sizes,each=n)
  dist <- dnorm(e,mu[2],sigma[4])
  dist <- dist/sum(dist)
  dist <- rep(dist,times=n)
  Z.res[,i] <- dist*sizes
}
run.b <- Z.res ## - these lines can be ignored. Was saving outputs for SI figure.
# calculate moments
moment.res <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
for (i in 1:n.gens){
  moment.res[,i] <- moments.fun(Z.res[,i],G)[1:4]
  moment.res.z[,i] <- moments.fun(Z.res[,i],Z)[1:4]
  moment.res.e[,i] <- moments.fun(Z.res[,i],E)[1:4]
}
# now plot the second figre
par(mfg=c(1,2))
plot(1:n.gens,moment.res[1,1:n.gens],type='l',xlab='',ylab='',ylim=c(0,30),axes=FALSE)#ylim=c(ranger[1],ranger[2])
box()
lines(1:n.gens,moment.res.z[1,1:n.gens],col='red')
lines(1:n.gens,moment.res.e[1,1:n.gens],col='blue')
axis(2,labels=TRUE)
axis(1,labels=TRUE)
mtext(side=2,line=2.5,'Mean')
text(8,28,'(b)',cex=1.5)

par(mfg=c(2,2))
plot(1:n.gens, (apply(Z.res,2,sum)[1:n.gens]),type='l',xlab='',ylab='',axes=FALSE)
mtext(side=2,line=2.5,'N')
axis(2,labels=TRUE)
axis(1,labels=TRUE)
box()
text(8,180,'(e)',cex=1.5)
# now ove onto the next simulation
####################################
Z.res.fit <- Z.res <- array(NA,c(n*n,n.gens))
# distribution of breeding values, environmental component and phenotype at time 1
Z.res[,1] <- dmvnorm(Z.cols,mean=mu,sigma=sigma)
Z.res[,1] <- Z.res[,1]/sum(Z.res[,1])
# parameters for fitness function -- same fitness funciton params but now dd in inheritance of E
int <- 0.3
slp.z1 <- 0.1
slp.N <- -0.01
# Run the simulation 
for (i in 2:n.gens){
  popN <- sum(Z.res[,i-1]) # calculate pop size
  w <- fitness.fun(G+E,int,slp.z1,slp.N,popN) # estimate fitness
  Z.res.fit[,i-1] <- w*Z.res[,i-1] # fitness in next time step
  sizes <- tapply(Z.res.fit[,i-1],G.count,sum) # fitness of G
  Es <- tapply(Z.res.fit[,i-1],E.count,sum) # Fitness of each E
  sizes <- rep(sizes,each=n)
  #    dist <- dnorm(e,mu[2],sigma[4])
  mat <- t(outer(e,e,G.fun,19,0,-0.065,1,popN)) #16,0.7,-0.07 # 20,0,-0.07
  mat <- mat/(matrix(apply(mat,2,sum),n,n,byrow=TRUE))
  dist <- mat%*%Es
  dist <- dist/sum(dist)
  dist <- rep(dist,times=n)
  Z.res[,i] <- dist*sizes
}
run.c <- Z.res
# calculate the moments
moment.res <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
for (i in 1:n.gens){
  moment.res[,i] <- moments.fun(Z.res[,i],G)[1:4]
  moment.res.z[,i] <- moments.fun(Z.res[,i],Z)[1:4]
  moment.res.e[,i] <- moments.fun(Z.res[,i],E)[1:4]
}
# And draw the graphs for the model with adaptive plasticity
par(mfg=c(1,3))
plot(1:n.gens,moment.res[1,1:n.gens],type='l',xlab='',ylab='',ylim=c(0,30),axes=FALSE) #ylim=c(ranger[1],ranger[2]),
box()
lines(1:n.gens,moment.res.z[1,1:n.gens],col='red')
lines(1:n.gens,moment.res.e[1,1:n.gens],col='blue')
axis(side=2,labels=TRUE)
mtext(side=2,line=2.5,'Mean')
axis(side=1,labels=TRUE)
text(8,28,'(c)',cex=1.5)
par(mfg=c(2,3))
plot(1:n.gens, (apply(Z.res,2,sum)[1:n.gens]),type='l',xlab='',ylab='',axes=FALSE)
box()
axis(side=2,labels=TRUE)
mtext(side=2,line=2.5,'N')
axis(side=1,labels=TRUE)
text(95,10,'(f)',cex=1.5)


#####
# the next block of code runs models with no plasticity, adaptive plasticity and non-adaptive plasticity
n.gens <- 100

#########################################

Z.res.fit <- Z.res <- array(NA,c(n*n,n.gens))
# distribution of breeding values, environmental component and phenotype at time 1
Z.res[,1] <- dmvnorm(Z.cols,mean=mu,sigma=sigma)
Z.res[,1] <- Z.res[,1]/sum(Z.res[,1])*116

# parameters for fitness function
int <- 2.16
slp.z1 <- 0
slp.N <- -0.01
# run the first simulation.
for (i in 2:n.gens){
  if (i==20) slp.z1 <- 0.1
  if (i==20){
    ggg <- moments.fun(Z.res[,i-1],Z)
    int <- 0.3
    temp.temp <- Z.res[,i-1]
    Z.res[,i-1] <- temp1 <- ifelse(Z>rep(ggg[1],n*n),Z.res[,i-1]*0.01,Z.res[,i-1])
  }
  popN <- sum(Z.res[,i-1]) # calculate pop size
  w <- fitness.fun(G+E,int,slp.z1,slp.N,popN) # estimate fitness
  Z.res.fit[,i-1] <- w*Z.res[,i-1] # fitness in next time step
  for (j in 1:n){
    dist <- Z.res.fit[G==g[j],i-1]
    mat <- t(outer(e,e,G.fun,12,0,0,1,popN)) #16,0.7,-0.07 # 20,0,-0.07
    mat <- mat/(matrix(apply(mat,2,sum),n,n,byrow=TRUE))
    dist <- mat%*%dist
    Z.res[G==g[j],i] <- dist
  }
  cat(paste(i,' '))# this whole section could do with being optimised, but given I only need to run the code once I really coulnd't be arsed.
}
run.d <- Z.res # ignore
# calculate moments of results distribution
moment.res <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
for (i in 1:n.gens){
  moment.res[,i] <- moments.fun(Z.res[,i],G)[1:4]
  moment.res.z[,i] <- moments.fun(Z.res[,i],Z)[1:4]
  moment.res.e[,i] <- moments.fun(Z.res[,i],E)[1:4]
}
ranger <- range(moment.res[1,],moment.res.z[1,],moment.res.e[1,1:n.gens])
# plot results from model with no plasticity
par(mfg=c(3,1))
plot(1:50,moment.res[1,1:50],type='l',xlab='',ylab='',ylim=c(5,25),axes=FALSE)
box()
axis(2,labels=TRUE)
mtext(side=2,line=2.5,'Means')
axis(1,labels=TRUE)
lines(1:n.gens,moment.res.z[1,],col='red')
lines(1:n.gens,moment.res.e[1,],col='blue')
temp <- c(moment.res[2,],moment.res.z[2,])
text(5,23,'(g)',cex=1.5)
par(mfg=c(4,1))
plot(1:50, (apply(Z.res,2,sum)[1:50]),type='l',xlab='',ylab='',axes=FALSE,ylim=c(60,140))
box()
mtext(side=1,line=2.5,'Time')
axis(1,labels=TRUE)
axis(2,labels=TRUE)
mtext(side=2,line=2.5,'N')
text(5,130,'(j)',cex=1.5)

#######################
# Now repeat the whole exercise for a model with adaptive plasticity. 
aaa <- 0.03

Z.res.fit <- Z.res <- array(NA,c(n*n,n.gens))

# distribution of breeding values, environmental component and phenotype at time 1
Z.res[,1] <- dmvnorm(Z.cols,mean=mu,sigma=sigma)
Z.res[,1] <- Z.res[,1]/sum(Z.res[,1])*116

# parameters for fitness function
int <- 2.16
slp.z1 <- 0
slp.N <- -0.01

for (i in 2:n.gens){
  if (i==20) slp.z1 <- 0.1
  if (i==20){
    ggg <- moments.fun(Z.res[,i-1],Z)
    int <- 0.3
    Z.res[,i-1] <- temp1 <- ifelse(Z>rep(ggg[1],n*n),Z.res[,i-1]*0.01,Z.res[,i-1])
  }
  popN <- sum(Z.res[,i-1]) # calculate pop size
  w <- fitness.fun(G+E,int,slp.z1,slp.N,popN) # estimate fitness
  Z.res.fit[,i-1] <- w*Z.res[,i-1] # fitness in next time step
  for (j in 1:n){
    dist <- Z.res.fit[G==g[j],i-1]
    mat <- t(outer(e,e,G.fun,12+aaa*116,0,-aaa,1,popN)) #16,0.7,-0.07 # 20,0,-0.07
    mat <- mat/(matrix(apply(mat,2,sum),n,n,byrow=TRUE))
    dist <- mat%*%dist
    Z.res[G==g[j],i] <- dist
  }
  cat(paste(i,' ')) # same deal, because I only run this once I am being lazy.  I want functionality, not elegance
}
run.e <- Z.res
moment.res <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
for (i in 1:n.gens){
  moment.res[,i] <- moments.fun(Z.res[,i],G)[1:4]
  moment.res.z[,i] <- moments.fun(Z.res[,i],Z)[1:4]
  moment.res.e[,i] <- moments.fun(Z.res[,i],E)[1:4]
}
ranger <- range(moment.res[1,1:50],moment.res.z[1,1:50],moment.res.e[1,1:50])
par(mfg=c(3,2))
plot(1:50,moment.res[1,1:50],type='l',xlab='',ylab='',ylim=c(5,25),axes=FALSE)
box()
axis(side=1,labels=TRUE)
axis(side=2,labels=TRUE)
mtext(side=2,line=2.5,'Means')
lines(1:50,moment.res.z[1,1:50],col='red')
lines(1:50,moment.res.e[1,1:50],col='blue')
temp <- c(moment.res[2,1:50],moment.res.z[2,])
text(5,23,'(h)',cex=1.5)
par(mfg=c(4,2))
plot(1:50, (apply(Z.res,2,sum)[1:50]),type='l',xlab='',ylab='',axes=FALSE,ylim=c(60,140))
axis(1,labels=TRUE)
axis(2,labels=TRUE)
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'N')
text(5,130,'(k)',cex=1.5)
box()

#############
# and now repeat the whole thing over again for the final model with non-adaptive plasticity.
# again -- being lazy as this code will only ever be run once.
Z.res.fit <- Z.res <- array(NA,c(n*n,n.gens))

# distribution of breeding values, environmental component and phenotype at time 1
Z.res[,1] <- dmvnorm(Z.cols,mean=mu,sigma=sigma)
Z.res[,1] <- Z.res[,1]/sum(Z.res[,1])*116

# parameters for fitness function
int <- 2.16
slp.z1 <- 0
slp.N <- -0.01

for (i in 2:n.gens){
  if (i==20) slp.z1 <- 0.1
  if (i==20){
    #    slp.z1 <- 0 else slp.z1 <- 0.01 # this line imposes selection at time point 150
    ggg <- moments.fun(Z.res[,i-1],Z)
    int <- 0.3
    Z.res[,i-1] <- temp1 <- ifelse(Z>rep(ggg[1],n*n),Z.res[,i-1]*0.01,Z.res[,i-1])
  }
  popN <- sum(Z.res[,i-1]) # calculate pop size
  w <- fitness.fun(G+E,int,slp.z1,slp.N,popN) # estimate fitness
  Z.res.fit[,i-1] <- w*Z.res[,i-1] # fitness in next time step
  for (j in 1:n){
    dist <- Z.res.fit[G==g[j],i-1]
    mat <- t(outer(e,e,G.fun,12-aaa*116,0,aaa,1,popN)) #16,0.7,-0.07 # 20,0,-0.07
    mat <- mat/(matrix(apply(mat,2,sum),n,n,byrow=TRUE))
    dist <- mat%*%dist
    Z.res[G==g[j],i] <- dist
  }
  cat(paste(i,' '))
}
run.f <- Z.res
moment.res <- moment.res.z <- moment.res.e <- array(NA,c(4,n.gens))
for (i in 1:n.gens){
  moment.res[,i] <- moments.fun(Z.res[,i],G)[1:4]
  moment.res.z[,i] <- moments.fun(Z.res[,i],Z)[1:4]
  moment.res.e[,i] <- moments.fun(Z.res[,i],E)[1:4]
}
ranger <- range(moment.res[1,n.gens],moment.res.z[1,n.gens],moment.res.e[1,n.gens])
par(mfg=c(3,3))
plot(1:50,moment.res[1,1:50],type='l',xlab='',ylab='',ylim=c(5,25),axes=FALSE)
box()
leg.txt <- c('Phenotype','Breeding value','Env component')
leg.lty <- c(1,1,1)
leg.col <- c('red','black','blue')
legend('topleft',legend=leg.txt,lty=leg.lty,col=leg.col,bty='n',cex=0.9)
mtext(side=2,line=2.5,'Means')
axis(1,labels=TRUE)
axis(2,labels=TRUE)

lines(1:50,moment.res.z[1,1:50],col='red')
lines(1:50,moment.res.e[1,1:50],col='blue')
temp <- c(moment.res[2,1:50],moment.res.z[2,])
text(45,23,'(i)',cex=1.5)
par(mfg=c(4,3))
plot(1:50, (apply(Z.res,2,sum)[1:50]),type='l',xlab='',ylab='',axes=FALSE,ylim=c(60,140))
axis(1,labels=TRUE)
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'N')
axis(2,labels=TRUE)
text(5,130,'(l)',cex=1.5)
box()



