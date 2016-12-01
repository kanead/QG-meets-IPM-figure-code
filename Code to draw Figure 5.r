# This code draws figure 5 in Coulson et al.
# It runs two simulations to examine how genetic covariances influence evolution (a) to (h)
# before examining how fitness functions would evolve for each character if the other
# character was not measured.

# clear memory
rm(list=ls())
# load key libraries
library(mvtnorm)
library(Matrix)
# close any open graphics
graphics.off()
# phenomenological inheritance function H(z'|z,t)  zz = z' for R code
G.fun <- function(z, zz, mu.a, mu.b, sigma.a) {
  mu.z <- mu.a+mu.b*z
  sigma.z2 <- sigma.a
  sigma.z <- sqrt(sigma.z2)
  temp1 <- sqrt(2*pi)*sigma.z
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  return(exp(-temp2)/temp1)
}
# calcaulte the moments of a distribution 
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
# open drawing canvas
quartz()
par(mfrow=c(3,4),mar=c(5,2,1,1),oma=c(3,3,1,1))

# set up global parameters and choose initial values
# number of classes in the matrix
n <- 100
# values of breeding value distribution
z <- seq(4,20,length.out=n)
x <- z
# create the combinations of x and y
z1 <- rep(z,each=n)
z2 <- rep(z,times=n)
h <- z2[2]-z2[1]

z <- cbind(z1,z2)
rm(z1,z2)
# number of generations for the simulation
n.gens <- 120
# define arrays into which to place results
A <- array(NA,c(n*n,n.gens))
# start with the initial distributions
mu.a1 <- c(7,15)
sigma.a1 <- matrix(c(1,-0.7,-0.7,1),2,2)
# distribution of breeding values, environmental component and phenotype at time 1
A[,1] <- dmvnorm(z,mean=mu.a1,sigma=sigma.a1)

# scale to sum to unity
A[,1] <- A[,1]/sum(A[,1])

# fitness function
# parameters for fitness function
int <- 2 #1.2
slp.z1 <- 0.15 
slp.z2 <- -0.13
fitness.fun <- function(x1,x2,int,slp1,slp2) int+slp1*x1+slp2*x2

# start the main loop here
for (gens in 2:n.gens){
  w <- fitness.fun(z[,1],z[,2],int,slp.z1,slp.z2)
  # breeding value fitness
  A[,gens] <- w*A[,gens-1]
}

# now to graphing
covar <- function(x,y,n) sum(x*y*n)/sum(n)-((sum(x*n)/sum(n))*(sum(y*n)/sum(n)))
means <- function(x,n) sum(x*n)/sum(n)
means.a1 <- means.a2 <- var.a1 <- var.a2 <- covar.a <- means.w <- rep(NA,n.gens)
for (i in 1:n.gens){
  var.a1[i] <- covar(z[,1],z[,1],A[,i])
  var.a2[i] <- covar(z[,2],z[,2],A[,i])
  covar.a[i] <- covar(z[,1],z[,2],A[,i])
  means.a1[i] <- means(z[,1],A[,i])
  means.a2[i] <- means(z[,2],A[,i])
  means.w[i] <- means(w,A[,i])
}
contour(x,x,matrix(A[,1]+A[,n.gens]/sum(A[,gens]),100,100),xlab='',ylab='',axes=FALSE)
axis(2,labels=TRUE)
axis(1,labels=TRUE)
box()
mtext(side=2,line=2.5,'Character 2')
mtext(side=1,line=2.5,'Character 1',col='red')
mtext(side=3,line=0.5,'Trait dynamics')
text(10,6,'Time 1',cex=1.5)
text(13,17,'Time 120',cex=1.5)
text(6,19,'(a)',cex=1.5)
lims <- range(c(means.a1,means.a2))
plot(1:n.gens,means.a1,type='l',ylim=c(lims[1],lims[2]),xlab='',ylab='',axes=FALSE)
lines(1:n.gens,means.a2,col='red')
axis(side=2,labels=TRUE)
box()
mtext(side=3,line=0.5,'Means')
text(20,14.85,'(b)',cex=1.5)
mtext(side=1,line=2.5,'Time')
axis(1,labels=TRUE)
lims.v <- range(c(var.a1,var.a2))
plot(1:n.gens,var.a1,xlab='',ylab='',type='l',ylim=c(lims.v[1],lims.v[2]),axes=FALSE)
lines(1:n.gens,var.a2,col='red')
mtext(side=3,line=0.5,'Variances')
mtext(side=1,line=2.5,'Time')
axis(1,labels=TRUE)
box()
axis(2,labels=TRUE)
text(20,0.972,'(c)',cex=1.5)
leg.txt <- c('Char 1','Char 2')
leg.col <- c('red','black')
leg.lty <- c(1,1)
legend('topright',legend=leg.txt,col=leg.col,lty=leg.lty,bty='n')

plot(1:n.gens,covar.a,type='l',xlab='',ylab='',axes=FALSE)
mtext(side=3,line=0.5,'Covariance')
box()
axis(2,labels=TRUE)
text(16,-0.375,'(d)',cex=1.5)
mtext(side=1,line=2.5,'Time')
axis(1,labels=TRUE)

hat.int.2 <- hat.slps.2 <- hat.int.1 <- hat.slps.1 <- rep(NA,n.gens)
for (i in 1:n.gens){
  hat.slps.1[i] <- slp.z1+covar.a[i]/var.a1[i]*slp.z2
  hat.int.1[i] <- means.w[i]-hat.slps.1[i]*means.a1[i]
  hat.slps.2[i] <- slp.z2+covar.a[i]/var.a2[i]*slp.z1
  hat.int.2[i] <- means.w[i]-hat.slps.2[i]*means.a2[i]
}

biv.m1.list <- list(hat.slps.1,hat.int.1,hat.slps.2,hat.int.2)
means.m1 <- list(means.a1,means.a2)


# now repeat the exercise but with a different covariance between characters
n <- 100
# values of breeding value distribution
z <- seq(4,20,length.out=n)
x <- z
# create the combinations of x and y
z1 <- rep(z,each=n)
z2 <- rep(z,times=n)
h <- z2[2]-z2[1]

z <- cbind(z1,z2)
rm(z1,z2)
# number of generations for the simulation
n.gens <- 120
# define arrays into which to place results
A <- array(NA,c(n*n,n.gens))
# start with the initial distributions
mu.a1 <- c(7,15)
sigma.a1 <- matrix(c(1,0.7,0.7,1),2,2)
# distribution of breeding values, environmental component and phenotype at time 1
A[,1] <- dmvnorm(z,mean=mu.a1,sigma=sigma.a1)

# scale to sum to unity
A[,1] <- A[,1]/sum(A[,1])

# fitness function
# parameters for fitness function
int <- 2 #1.2
slp.z1 <- 0.15 
slp.z2 <- -0.13
fitness.fun <- function(x1,x2,int,slp1,slp2) int+slp1*x1+slp2*x2

# start the main loop here
for (gens in 2:n.gens){
  w <- fitness.fun(z[,1],z[,2],int,slp.z1,slp.z2)
  # breeding value fitness
  A[,gens] <- w*A[,gens-1]
}

# now to graphing - again
covar <- function(x,y,n) sum(x*y*n)/sum(n)-((sum(x*n)/sum(n))*(sum(y*n)/sum(n)))
means <- function(x,n) sum(x*n)/sum(n)
means.a1 <- means.a2 <- var.a1 <- var.a2 <- covar.a <- rep(NA,n.gens)
for (i in 1:n.gens){
  var.a1[i] <- covar(z[,1],z[,1],A[,i])
  var.a2[i] <- covar(z[,2],z[,2],A[,i])
  covar.a[i] <- covar(z[,1],z[,2],A[,i])
  means.a1[i] <- means(z[,1],A[,i])
  means.a2[i] <- means(z[,2],A[,i])
  means.w[i] <- means(w,A[,i])
}

contour(x,x,matrix(A[,1]+A[,n.gens]/sum(A[,gens]),100,100),xlab='',ylab='',axes=FALSE)
mtext(side=2,line=2.5,'Character 2')
box()
axis(side=2,labels=TRUE)
text(8,6,'Time 1',cex=1.5)
text(13,14,'Time 120',cex=1.5)
text(6,19,'(e)',cex=1.5)
mtext(side=1,line=2.5,'Character 1',col='red')
axis(1,labels=TRUE)

plot(1:n.gens,means.a1,type='l',ylim=c(lims[1],lims[2]),xlab='',ylab='',axes=FALSE)
lines(1:n.gens,means.a2,col='red')
#mtext(side=1,line=2.5,'Time')
box()
axis(side=2,labels=TRUE)
text(108,14.8,'(f)',cex=1.5)
mtext(side=1,line=2.5,'Time')
axis(1,labels=TRUE)

plot(1:n.gens,var.a1,xlab='',ylab='',type='l',ylim=c(lims.v[1],lims.v[2]),axes=FALSE)
lines(1:n.gens,var.a2,col='red')
#mtext(side=1,line=2.5,'Time')
box()
axis(side=2,labels=TRUE)
text(15,0.67,'(g)',cex=1.5)
mtext(side=1,line=2.5,'Time')
axis(1,labels=TRUE)

plot(1:n.gens,covar.a,type='l',xlab='',ylab='',axes=FALSE)
box()
axis(side=2,labels=TRUE)
text(15,0.733,'(h)',cex=1.5)
mtext(side=1,line=2.5,'TIme')
axis(1,labels=TRUE)

# now move to estimating fitness functions in the absence of each character being measured
hat.int.2 <- hat.slps.2 <- hat.int.1 <- hat.slps.1 <- rep(NA,n.gens)
for (i in 1:n.gens){
  hat.slps.1[i] <- slp.z1+covar.a[i]/var.a1[i]*slp.z2
  hat.int.1[i] <- means.w[i]-hat.slps.1[i]*means.a1[i]
  hat.slps.2[i] <- slp.z2+covar.a[i]/var.a2[i]*slp.z1
  hat.int.2[i] <- means.w[i]-hat.slps.2[i]*means.a2[i]
}

biv.m2.list <- list(hat.slps.1,hat.int.1,hat.slps.2,hat.int.2)
means.m2 <- list(means.a1,means.a2)



hat.int.2 <- hat.slps.2 <- hat.int.1 <- hat.slps.1 <- rep(NA,n.gens)
for (i in 1:n.gens){
  hat.slps.1[i] <- slp.z1+covar.a[i]/var.a1[i]*slp.z2
  hat.int.1[i] <- means.w[i]-hat.slps.1[i]*means.a1[i]
  hat.slps.2[i] <- slp.z2+covar.a[i]/var.a2[i]*slp.z1
  hat.int.2[i] <- means.w[i]-hat.slps.2[i]*means.a2[i]
}

biv.m3.list <- list(hat.slps.1,hat.int.1,hat.slps.2,hat.int.2)
means.m3 <- list(means.a1,means.a2)

plot(x,biv.m1.list[[2]][1]+biv.m1.list[[1]][1]*x,type='l',xlab='',ylab='')
lines(x,biv.m1.list[[2]][n.gens]+biv.m1.list[[1]][n.gens]*x,lty=2)
lines(x,biv.m2.list[[2]][1]+biv.m2.list[[1]][1]*x,col='red')
lines(x,biv.m2.list[[2]][n.gens]+biv.m2.list[[1]][n.gens]*x,col='red',lty=2)
leg.text <- c('-cov','-cov')
leg.lty <- c(1,1)
leg.col <- c('black','red')
legend('topleft',legend=leg.text,lty=leg.lty,col=leg.col,bty='n')
mtext(side=3,line=0.5,'')
mtext(side=1,line=2.5,'Character 1')
mtext(side=2,line=2.5,'Fitness')
text(18.8,0.7,'(i)',cex=1.5)
mtext(side=3,line=0.5,'Fitness function')
text(18,1.45,'t=1')
text(18,2.5,'t=120')
plot(x,biv.m1.list[[4]][1]+biv.m1.list[[3]][1]*x,type='l',xlab='',ylab="")
lines(x,biv.m1.list[[4]][n.gens]+biv.m1.list[[3]][n.gens]*x,lty=2)
lines(x,biv.m2.list[[4]][1]+biv.m2.list[[3]][1]*x,col='red')
lines(x,biv.m2.list[[4]][n.gens]+biv.m2.list[[3]][n.gens]*x,col='red',lty=2)

text(5.3,0.2,'(j)',cex=1.5)
mtext(side=3,line=0.5,'Fitness function')
mtext(side=1,line=2.5,'Character 2')

sel.m1.z1 <- means.m1[[1]][-1]-means.m1[[1]][-length(means.m1[[1]])]
sel.m1.z2 <- means.m1[[2]][-1]-means.m1[[2]][-length(means.m1[[1]])]
sel.m2.z1 <- means.m2[[1]][-1]-means.m2[[1]][-length(means.m2[[1]])]
sel.m2.z2 <- means.m2[[2]][-1]-means.m2[[2]][-length(means.m2[[1]])]
sel.m3.z1 <- means.m3[[1]][-1]-means.m3[[1]][-length(means.m3[[1]])]
sel.m3.z2 <- means.m3[[2]][-1]-means.m3[[2]][-length(means.m3[[1]])]
