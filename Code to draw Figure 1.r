# This file draws Figure 1 in Coulson et al. Modeling Adaptive and Non-adaptive Responses of Populations to Environmental Change
# clear memory and shut any graphics windows
# The code is not elegant. It is functional.  Could wrap the repeated graphics call into a function quite easily but life is sometimes too short.
rm(list=ls())
graphics.off()
# load appropriate libraries
library(mvtnorm)
library(Matrix)

# define a linear fitness function
w <- function(z,a,b) (a+b*z)

# a function to calculate moments of a distribution within an IPM. x is the density 
# and z is the mid point trait values.
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

#covariance function - x and y are values of the trait and n is the density of each combination
covar <- function(x,y,n) sum(x*y*n)/sum(n)-((sum(x*n)/sum(n))*(sum(y*n)/sum(n)))

# inheritance funciton - calculates the probability of transition from z to zz (X to X' in paper)
H.fun <- function(z, zz, mu.a, mu.b, sigma.a) {
  mu.z <- mu.a+mu.b*z
  sigma.z2 <- sigma.a
  sigma.z <- sqrt(sigma.z2)
  temp1 <- sqrt(2*pi)*sigma.z
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  return(exp(-temp2)/temp1)
}
# fitness functions
surv.fun <- function(z1,int,slp) 1/(1+exp(-(int+slp*z1)))
repr.fun <- function(z1,int,slp) exp(int+slp*z1)

# number of bins -- how many bins to discretise the matrix approximation of the IPM into
n <- 500

# number of iterations -- how long to run the simulation for
n.gens <- 30

# midpoints -- these are the values of the character
z <- seq(0,30,length.out=n)

# results vector -- where to store the densities of the distributions
res.p <- array(NA,c(n,n.gens))

# starting condition -- start with a normal distribution scaled to sum to unity
res.p[,1] <- dnorm(z,8,1)
res.p[,1] <- res.p[,1]/sum(res.p[,1])

# calculate fitness of each value of the character, z
fit <- w(z,0.1,0.2)

# iterate forward the population
for (i in 2:n.gens){
  res.p[,i] <- fit*res.p[,i-1]
}

# calculate the population growth rate R0 -- first calculate population sizes, then growth rate
Ns <- apply(res.p,2,sum)
R0 <- Ns[-1]/Ns[-length(Ns)]

# calculaate moments of the distribution at each generation
mom.res.sel <- array(NA,c(n.gens,4))
for (i in 1:n.gens) mom.res.sel[i,] <- moments.fun(res.p[,i],z)[1:4]


# we now repeat the process imposing a genetic constraint
# define an array in which to store results
res.p.c <- array(NA,c(n,n.gens))

# starting condition -- no values of z greater than 11.5
res.p.c[,1] <- dnorm(z,8,1)
res.p.c[,1] <- ifelse(z>11.5,0,res.p.c[,1])
res.p.c[,1] <- res.p.c[,1]/sum(res.p.c[,1])

# iteracte forward the popuulation as above
for (i in 2:n.gens){
  res.p.c[,i] <- fit*res.p.c[,i-1]
}

# now calculate the population size, population growth rate and moments of the distribution
Ns.c <- apply(res.p.c,2,sum)
R0.c <- Ns.c[-1]/Ns.c[-length(Ns.c)]

mom.res.sel.c  <- array(NA,c(n.gens,4))
for (i in 1:n.gens) mom.res.sel.c[i,] <- moments.fun(res.p.c[,i],z)[1:4]

# define a place to store the graphics
quartz() #png('~/dropbox/Inheritance paper/Paper 1 in series/Figure3.png')

# This first block of code compares dynamics for models with and without a genetic constraint
par(mfrow=c(3,4),mar=c(4,2,2,1),oma=c(0,3,0,0))

plot(1:(n.gens-1),R0,xlab='',ylab='',type='l',xlim=c(1,n.gens-2),axes=FALSE)
lines(1:(n.gens-1),R0.c,col='red')
box()
axis(2,labels=TRUE)
axis(1,labels=TRUE)
mtext(side=1,line=2.5,'Time')
mtext(side=3,line=0.5,'Population growth')
text(4,2.175,'(a)',cex=1.5)

plot(1:n.gens,mom.res.sel[,1],xlab='',ylab='',type='l',xlim=c(1,n.gens-2),axes=FALSE)
lines(1:n.gens,mom.res.sel.c[,1],col='red')
box()
axis(2,labels=TRUE)
mtext(side=3,line=0.5,'Mean')
axis(1,labels=TRUE)
mtext(side=1,line=2.5,'Time')
text(4,10.45,'(b)',cex=1.5)

temp <- range(c(mom.res.sel[,2],mom.res.sel.c[,2]))
plot(1:n.gens,mom.res.sel[,2],xlab='',ylab='',type='l',xlim=c(1,n.gens-2),ylim=c(temp[1],temp[2]),axes=FALSE)
lines(1:n.gens,mom.res.sel.c[,2],col='red')
box()
axis(2,labels=TRUE)
mtext(side=3,line=0.5,'Variance')
text(25,0.965,'(c)',cex=1.5)
axis(1,labels=TRUE)
mtext(side=1,line=2.5,'Time')

temp <- range(c(mom.res.sel[,3],mom.res.sel.c[,3]))
plot(1:n.gens,mom.res.sel[,3],xlab='',ylab='',type='l',xlim=c(1,n.gens-2),ylim=c(temp[1],temp[2]),axes=FALSE)
lines(1:n.gens,mom.res.sel.c[,3],col='red')
box()
axis(2,labels=TRUE)
mtext(side=3,line=0.5,'Skew')
text(4,-0.53,'(d)',cex=1.5)
axis(1,labels=TRUE)
mtext(side=1,line=2.5,'Time')

# We now run a model for a bivariate character distribution
# set up global parameters and choose initial values
n <- 100 # number of classes in the matrix
z <- seq(4,20,length.out=n) # values of breeding value distribution of bivariate distribution
x <- z
z2 <- rep(z,each=n); z1 <- rep(z,times=n); z <- cbind(z1,z2); rm(z1,z2)
n.gens <- 300 # number of generations for the simulation

# define arrays into which to place results
age.1 <- age.2 <- array(NA,c(n*n,n.gens))
# start with the initial distributions
mu.a1 <- c(7,10)
sigma.a1 <- matrix(c(1,0.7,0.7,0.8),2,2)
# distribution of breeding values, environmental component and phenotype at time 1
age.2[,1] <- age.1[,1] <- dmvnorm(z,mean=mu.a1,sigma=sigma.a1)

# scale to sum to unity
age.1[,1] <- age.1[,1]/sum(age.1[,1]); age.2[,1] <- age.2[,1]/sum(age.2[,1])

# now iterate the bivariate distribution forwards in time.
for (i in 1:(n.gens-1)){
  # selection operates on the trait at age 1
  age.2[,i+1] <- surv.fun(z[,1],0.1,0.06)*age.1[,i]
  # recruitment of age 2 individuals to construct age 1 individuals
  age.1[,i+1] <- repr.fun(z[,2],0.01,0.05)*age.2[,i]
}

# calculate covariances and moments of each distribution 
num <- covar(z[,1],z[,2],age.2[,1])
z1 <- tapply(age.2[,1],z[,1],sum)
z2 <- tapply(age.2[,1],z[,2],sum)
moms.1 <- moments.fun(z1,x)
moms.2 <- moments.fun(z2,x)
slp <- num/moms.1[2]
int <- moms.2[1]-moms.1[1]*slp
int.v <- moms.1[2]*(1-slp^2)
# now use these values to construct the phenomenological development function from age a to age a+1 
mat <- t(outer(x,x,H.fun,int,slp,int.v))
mat <- mat/matrix(apply(mat,2,sum),nrow=n,ncol=n,byrow=TRUE)

# now plot the distributions overlaid with the phenomenological maps figure (e).
x.v <- which(x>5 & x<10)
y.v <- which(x>7.5 & x<12.5)
mat <- t(mat)
image(x[x.v],x[y.v],mat[x.v,y.v],xlim=c(5,20),ylim=c(5,20),xlab='',ylab='',axes=FALSE)
box()
axis(2,labels=TRUE)
mtext(side=2,line=2.5,'Trait age 2')
axis(1,labels=TRUE)
mtext(side=1,line=2.5,'Trait age 1')
text(18,7,'(e)',cex=1.5)

num <- covar(z[,1],z[,2],age.2[,200])
z1 <- tapply(age.2[,200],z[,1],sum)
z2 <- tapply(age.2[,200],z[,2],sum)
moms.1 <- moments.fun(z1,x)
moms.2 <- moments.fun(z2,x)
slp <- num/moms.1[2]
int <- moms.2[1]-moms.1[1]*slp
int.v <- moms.1[2]*(1-slp^2)
mat <- t(outer(x,x,H.fun,int,slp,int.v))
mat <- mat/matrix(apply(mat,2,sum),nrow=n,ncol=n,byrow=TRUE)
x.v <- which(x>9.5 & x<14.5)
y.v <- which(x>12.5 & x<17.5)
mat.p <- t(mat)
par(new=TRUE)
image(x[x.v],x[y.v],mat.p[x.v,y.v],xlim=c(5,20),ylim=c(5,20),xlab='',ylab='',axes=FALSE)

par(new=TRUE)
temp <- matrix(age.1[,1]/sum(age.1[,1]),nrow=n,ncol=n)+matrix(age.1[,200]/sum(age.1[,200]),nrow=n,ncol=n)
contour(x,x,temp,xlim=c(5,20),ylim=c(5,20),xaxs='i',yaxs='i',axes=FALSE)
text(7.5,6.5,'Time 1')
text(11.7,18.5,'Time 200')


# Repeat whole exercise for a different initial covariance - figure (f)
# start with the initial distributions
mu.a1 <- c(7,10)
sigma.a1 <- matrix(c(1,-0.2,-0.2,0.8),2,2)
# distribution of breeding values, environmental component and phenotype at time 1
age.2[,1] <- age.1[,1] <- dmvnorm(z,mean=mu.a1,sigma=sigma.a1)

# scale to sum to unity
age.1[,1] <- age.1[,1]/sum(age.1[,1]); age.2[,1] <- age.2[,1]/sum(age.2[,1])

for (i in 1:(n.gens-1)){
  # selection operates on the trait at age 1
  age.2[,i+1] <- surv.fun(z[,1],0.1,0.06)*age.1[,i]
  age.1[,i+1] <- repr.fun(z[,2],0.01,0.05)*age.2[,i]
}

num <- covar(z[,1],z[,2],age.2[,1])
z1 <- tapply(age.2[,1],z[,1],sum)
z2 <- tapply(age.2[,1],z[,2],sum)
moms.1 <- moments.fun(z1,x)
moms.2 <- moments.fun(z2,x)
slp <- num/moms.1[2]
int <- moms.2[1]-moms.1[1]*slp
int.v <- moms.1[2]*(1-slp^2)
mat <- t(outer(x,x,H.fun,int,slp,int.v))
mat <- mat/matrix(apply(mat,2,sum),nrow=n,ncol=n,byrow=TRUE)

x.v <- which(x>5 & x<10)
y.v <- which(x>7.5 & x<12.5)
mat <- t(mat)
image(x[x.v],x[y.v],mat[x.v,y.v],xlim=c(5,20),ylim=c(5,20),xlab='',ylab='',axes=FALSE)
text(18,7,'(f)',cex=1.5)
mtext(side=1,line=2.5,'Trait age 1')
axis(1,labels=TRUE)
axis(2,labels=TRUE)

num <- covar(z[,1],z[,2],age.2[,300])
z1 <- tapply(age.2[,300],z[,1],sum)
z2 <- tapply(age.2[,300],z[,2],sum)
moms.1 <- moments.fun(z1,x)
moms.2 <- moments.fun(z2,x)
slp <- num/moms.1[2]
int <- moms.2[1]-moms.1[1]*slp
int.v <- moms.1[2]*(1-slp^2)
mat <- t(outer(x,x,H.fun,int,slp,int.v))
mat <- mat/matrix(apply(mat,2,sum),nrow=n,ncol=n,byrow=TRUE)
x.v <- which(x>6 & x<11)
y.v <- which(x>13 & x<18)
mat.p <- t(mat)
par(new=TRUE)
image(x[x.v],x[y.v],mat.p[x.v,y.v],xlim=c(5,20),ylim=c(5,20),xlab='',ylab='',axes=FALSE)

par(new=TRUE)
temp <- matrix(age.1[,1]/sum(age.1[,1]),nrow=n,ncol=n)+matrix(age.1[,n.gens]/sum(age.1[,n.gens]),nrow=n,ncol=n)
contour(x,x,temp,xlim=c(5,20),ylim=c(5,20),xaxs='i',yaxs='i',axes=FALSE)
text(7.5,6.5,'Time 1')
text(8.5,18.75,'Time 300')
box()




# repeat exercise for a different set of values. Figure (g)
# start with the initial distributions
mu.a1 <- c(7.5,16)
sigma.a1 <- matrix(c(1,-0.1,-0.1,0.8),2,2)
# distribution of breeding values, environmental component and phenotype at time 1
age.2[,1] <- age.1[,1] <- dmvnorm(z,mean=mu.a1,sigma=sigma.a1)

# scale to sum to unity
age.1[,1] <- age.1[,1]/sum(age.1[,1]); age.2[,1] <- age.2[,1]/sum(age.2[,1])

for (i in 1:(n.gens-1)){
  # selection operates on the trait at age 1
  age.2[,i+1] <- surv.fun(z[,1],0.1,0.03)*age.1[,i]
  age.1[,i+1] <- repr.fun(z[,2],0.01,-0.075)*age.2[,i]
}

num <- covar(z[,1],z[,2],age.2[,1])
z1 <- tapply(age.2[,1],z[,1],sum)
z2 <- tapply(age.2[,1],z[,2],sum)
moms.1 <- moments.fun(z1,x)
moms.2 <- moments.fun(z2,x)
slp <- num/moms.1[2]
int <- moms.2[1]-moms.1[1]*slp
int.v <- moms.1[2]*(1-slp^2)
mat <- t(outer(x,x,H.fun,int,slp,int.v))
mat <- mat/matrix(apply(mat,2,sum),nrow=n,ncol=n,byrow=TRUE)

x.v <- which(x>5 & x<10)
y.v <- which(x>13.5 & x<18.5)
mat <- t(mat)
image(x[x.v],x[y.v],mat[x.v,y.v],xlim=c(5,20),ylim=c(5,20),xlab='',ylab='')
text(18,7,'(g)',cex=1.5)
num <- covar(z[,1],z[,2],age.2[,300])
z1 <- tapply(age.2[,300],z[,1],sum)
z2 <- tapply(age.2[,300],z[,2],sum)
moms.1 <- moments.fun(z1,x)
moms.2 <- moments.fun(z2,x)
slp <- num/moms.1[2]
int <- moms.2[1]-moms.1[1]*slp
int.v <- moms.1[2]*(1-slp^2)
mat <- t(outer(x,x,H.fun,int,slp,int.v))
mat <- mat/matrix(apply(mat,2,sum),nrow=n,ncol=n,byrow=TRUE)
x.v <- which(x>8 & x<13)
y.v <- which(x>4.5 & x<9.5)
mat.p <- t(mat)
par(new=TRUE)
image(x[x.v],x[y.v],mat.p[x.v,y.v],xlim=c(5,20),ylim=c(5,20),xlab='',ylab='',axes=FALSE)

par(new=TRUE)
temp <- matrix(age.1[,1]/sum(age.1[,1]),nrow=n,ncol=n)+matrix(age.1[,n.gens]/sum(age.1[,n.gens]),nrow=n,ncol=n)
contour(x,x,temp,xlim=c(5,20),ylim=c(5,20),xaxs='i',yaxs='i',axes=FALSE)
text(8,19,'Time 1')
text(10.5,10.5,'Time 300')
mtext(side=1,line=2.5,'Trait age 1')


