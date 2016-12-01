# this code draws figure 3 in Coulson et al. 
# clear memory
# close graphics
# once again, the code is not elegant, but it is functional!
rm(list=ls())
graphics.off()
# open graphics window and set up drawing area
quartz()
par(mfrow=c(3,3),mar=c(4,4,1,1))

# function to calculate moments of a distribution
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

# Gaussian inheritance function H(z'|z,t)  zz = z' for R code
G.fun <- function(z, zz, mu.a, mu.b, sigma.a, n.b, N) {
  mu.z <- mu.a+mu.b*z+n.b*N
  sigma.z2 <- sigma.a
  sigma.z <- sqrt(sigma.z2)
  temp1 <- sqrt(2*pi)*sigma.z
  temp2 <- ((zz-mu.z)^2)/(2*sigma.z2)
  return(exp(-temp2)/temp1)
}

# fitness function
w <- function(z,a,b,c,N) a+b*z+c*N

# num bins
n <- 100

# num iterations
n.gens <- 500

# midpoints
z <- seq(0,30,length.out=n)

# results vector
res <- array(NA,c(n,n.gens))

# starting condition
res[,1] <- dnorm(z,8,1.1)
res[,1] <- res[,1]/sum(res[,1])

# now run first simulation for figure (a) -- examine consequences of a one of perturbation on recovery times
for (i in 2:n.gens){
  N <- sum(res[,i-1])
  fit <- w(z,0.2,0.1,-0.002,N)
  fit <- if (i==300) fit*0.4 else fit
  H <- t(outer(z,z,G.fun,5,0.5,1,-0.005,N))
  H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  res[,i] <- H%*%(fit*res[,i-1])
}

# plot figure (a)
Ns <- apply(res,2,sum)
plot(281:380,Ns[281:380],type='l',lwd=2,xlab='',ylab='',axes=FALSE,col='red')
axis(1,labels=seq(0,100,length.out=9),at=c(seq(280,380,length.out=9)))
axis(2,labels=TRUE)
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'N')
# extract carrying capacity
K <- Ns[n.gens]
text(290,32,'(a)',cex=1.5)
res <- array(NA,c(n,n.gens))

# now run the second simulation for figure (a)
# starting condition
res[,1] <- dnorm(z,8,1.1)
res[,1] <- res[,1]/sum(res[,1])

for (i in 2:n.gens){
  N <- sum(res[,i-1])
  fit <- w(z,0.2,0.1,-0.002,N)
  fit <- if (i==300) fit*0.4 else fit
  H <- t(outer(z,z,G.fun,5-0.005*K,0.5,1,0,N))
  H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  res[,i] <- H%*%(fit*res[,i-1])
}
Ns <- apply(res,2,sum)
lines(281:380,Ns[281:380],col='black',lwd=2)

box()

res <- array(NA,c(n,n.gens))

# now run third simulation for figure (a)
# starting condition
res[,1] <- dnorm(z,8,1.1)
res[,1] <- res[,1]/sum(res[,1])
for (i in 2:n.gens){
  N <- sum(res[,i-1])
  fit <- w(z,0.2,0.1,-0.002,N)
  fit <- if (i==300) fit*0.4 else fit
  H <- t(outer(z,z,G.fun,5-0.01*K,0.5,1,0.005,N))
  H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  res[,i] <- H%*%(fit*res[,i-1])
}
Ns <- apply(res,2,sum)
lines(281:380,Ns[281:380],col='blue',lwd=2)


# Now move onto figure (b)  This time have a permanent perturbation

# results vector
res <- array(NA,c(n,n.gens))

# first simulation for figure (b)
# starting condition
res[,1] <- dnorm(z,8,1.1)
res[,1] <- res[,1]/sum(res[,1])

for (i in 2:n.gens){
  N <- sum(res[,i-1])
  if (i>300) a <- 0.15 else a <- 0.2
  fit <- w(z,a,0.1,-0.002,N)
  H <- t(outer(z,z,G.fun,5,0.5,1,-0.005,N))
  H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  res[,i] <- H%*%(fit*res[,i-1])
}
res.1 <- res
Ns <- apply(res,2,sum)
K <- Ns[299]
plot(281:380,Ns[281:380],type='l',lwd=2,xlab='',ylab='',axes=FALSE,ylim=c(20,K),col='red')
axis(1,labels=seq(0,100,length.out=9),at=c(seq(280,380,length.out=9)))
axis(2,labels=TRUE)
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'N')
text(290,24,'(b)',cex=1.5)

res <- array(NA,c(n,n.gens))
# second simulation for figure (b)

# starting condition
res[,1] <- dnorm(z,8,1.1)
res[,1] <- res[,1]/sum(res[,1])

for (i in 2:n.gens){
  N <- sum(res[,i-1])
  if (i>300) a <- 0.15 else a <- 0.2
  fit <- w(z,a,0.1,-0.002,N)
  H <- t(outer(z,z,G.fun,5-0.005*K,0.5,1,0,N))
  H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  res[,i] <- H%*%(fit*res[,i-1])
}
Ns <- apply(res,2,sum)
lines(281:380,Ns[281:380],lwd=2)
res.2 <- res
box()

res <- array(NA,c(n,n.gens))

# Now third simulation for figure (c)
# starting condition
res[,1] <- dnorm(z,8,1.1)
res[,1] <- res[,1]/sum(res[,1])
for (i in 2:n.gens){
  N <- sum(res[,i-1])
  if (i>300) a <- 0.15 else a <- 0.2
  fit <- w(z,a,0.1,-0.002,N)
  H <- t(outer(z,z,G.fun,5-0.01*K,0.5,1,0.005,N))
  H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  res[,i] <- H%*%(fit*res[,i-1])
}
Ns <- apply(res,2,sum)
lines(281:380,Ns[281:380],col='blue',lwd=2)
res.3 <- res


# finally examine dynamics of mean phenotype for the three simulations in figure (b) and plot the means in Figure (c)
means.ph1 <- means.ph2 <- means.ph3 <- rep(NA,n.gens)
for (i in 1:500){
  means.ph1[i] <- moments.fun(res.1[,i],z)[1]
  means.ph2[i] <- moments.fun(res.2[,i],z)[1]
  means.ph3[i] <- moments.fun(res.3[,i],z)[1]
}
temp <- c(means.ph1[100:500],means.ph2[100:500],means.ph3[100:500])
plot(281:380,means.ph1[281:380],col='red',lwd=2,ylim=c(min(temp),max(temp)),type='l',axes=FALSE,xlab='',ylab='')
lines(281:380,means.ph2[281:380],col='black',lwd=2)
lines(281:380,means.ph3[281:380],col='blue',lwd=2)
axis(1,labels=seq(0,100,length.out=9),at=c(seq(280,380,length.out=9)))
axis(2,labels=TRUE)
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'Means')
box()
text(290,9.55,'(c)',cex=1.5)
