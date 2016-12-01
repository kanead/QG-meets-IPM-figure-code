# clear memory and turn off graphics
rm(list=ls())
graphics.off()

# open up drawing canvas (Quartz only works on a mac)
quartz()
par(mfrow=c(3,3),mar=c(4,4,1,1))

# specify number of generations to run simulation for
n.gens <- 10
# set up results vector
res <- frequencies <- array(NA,c(n.gens,3))
# fitness parameters -- additive effects of genotype
w <- c(0.85,1,1.15)
# initial conditions -- numbers of each gneotype
res[1,] <- c(10,10,10)
# run simulation of genotype dynamics -- note that inheritance is frequency dependent
for (i in 2:n.gens){
  sel <- res[i-1,]*w
  freq <- sel/sum(sel)
  prop.AA <- freq[1]^2 + freq[1]*freq[2] + 0.25*freq[2]^2
  prop.AB <- 0.5*freq[2]^2 + 2*freq[1]*freq[3] + freq[2]*freq[1] + freq[2]*freq[3]
  prop.BB <- freq[3]^2 + 0.25*freq[2]^2 + freq[2]*freq[3]
  gfreq <- c(prop.AA,prop.AB,prop.BB)
  res[i,] <- gfreq*sum(sel)
  frequencies[i,] <- gfreq
}

# values to plot.
new.cohort.A <- frequencies[n.gens,1]+0.5*frequencies[n.gens,2]
new.cohort.B <- frequencies[n.gens,3]+0.5*frequencies[n.gens,2]
after.sel.A <- freq[1]+0.5*freq[2]
after.sel.B <- freq[3]+0.5*freq[2]

# write out inheritance matrix for genotypes -- parents to offspring
mat <- matrix(NA,3,3)
mat[1,1] <- freq[1]^2 + 0.5*freq[1]*freq[2]
mat[2,1] <- 0.5*freq[1]*freq[2] + freq[1]*freq[3]
mat[3,1] <- 0

mat[1,2] <- 0.25*freq[2]^2 + 0.5*freq[1]*freq[2]
mat[2,2] <- 0.5*freq[1]*freq[2] + 0.5*freq[2]^2 + 0.5*freq[3]*freq[2]
mat[3,2] <- 0.25*freq[2]^2 + 0.5*freq[3]*freq[2]

mat[1,3] <- 0
mat[2,3] <- 0.5*freq[2]*freq[3] + freq[1]*freq[3]
mat[3,3] <- 0.5*freq[3]*freq[2] + freq[3]^2

# Plot Figure (e)
counts <- t(cbind(frequencies[n.gens-1,],freq,frequencies[n.gens,]))
plot(0,type='n',axes=FALSE,ann=FALSE)
par(mfg=c(2,2))
barplot(counts,beside=TRUE,col=c('red','black','blue'),xlab='',ylab='',ylim=c(0,1))
axis(1,labels=c('AA','AB','BB'),at=c(2.5,6.5,10.5))
box()
text(11,0.93,'(e)',cex=1.5)
mtext(side=1,line=2.5,'Genotype')
mtext(side=2,line=2.5,'Frequency')
leg.txt <- c('Before selection','After selection','After inheritance')
leg.col <- c('red','black','blue')
leg.lty <- c(1,1,1)
legend('topleft',col=leg.col,legend=leg.txt,lty=leg.lty,bty='n')

# now repeat process for the case when G is not additive with respective to g
res <- frequencies <- array(NA,c(100,3))
w <- c(0.65,1,0.35)
res[1,] <- c(10,10,10)

for (i in 2:100){
  sel <- res[i-1,]*w
  freq <- sel/sum(sel)
  prop.AA <- freq[1]^2 + freq[1]*freq[2] + 0.25*freq[2]^2
  prop.AB <- 0.5*freq[2]^2 + 2*freq[1]*freq[3] + freq[2]*freq[1] + freq[2]*freq[3]
  prop.BB <- freq[3]^2 + 0.25*freq[2]^2 + freq[2]*freq[3]
  gfreq <- c(prop.AA,prop.AB,prop.BB)
  res[i,] <- gfreq*sum(sel)
  frequencies[i,] <- gfreq
}

new.cohort.A <- frequencies[100,1]+0.5*frequencies[100,2]
new.cohort.B <- frequencies[100,3]+0.5*frequencies[100,2]
after.sel.A <- freq[1]+0.5*freq[2]
after.sel.B <- freq[3]+0.5*freq[2]

mat <- matrix(NA,3,3)
mat[1,1] <- freq[1]^2 + 0.5*freq[1]*freq[2]
mat[2,1] <- 0.5*freq[1]*freq[2] + freq[1]*freq[3]
mat[3,1] <- 0

mat[1,2] <- 0.25*freq[2]^2 + 0.5*freq[1]*freq[2]
mat[2,2] <- 0.5*freq[1]*freq[2] + 0.5*freq[2]^2 + 0.5*freq[3]*freq[2]
mat[3,2] <- 0.25*freq[2]^2 + 0.5*freq[3]*freq[2]

mat[1,3] <- 0
mat[2,3] <- 0.5*freq[2]*freq[3] + freq[1]*freq[3]
mat[3,3] <- 0.5*freq[3]*freq[2] + freq[3]^2

# now plot figure f -- non additive map from g to G
counts <- t(cbind(frequencies[100-1,],freq,frequencies[100,]))
par(mfg=c(2,3))
barplot(counts,beside=TRUE,col=c('red','black','blue'),xlab='',ylab='',ylim=c(0,0.6))
axis(1,labels=c('AA','AB','BB'),at=c(2.5,6.5,10.5))
box()
text(11,0.55,'(f)',cex=1.5)
mtext(side=1,line=2.5,'Genotype')
mtext(side=2,line=2.5,'Frequency')


# The next bit of code rather inelegantly projects
# a clonal model, a sexual model and our Gaussian approximation

# create a character of length 100
z <- seq(0,10,length.out=100)
# generate a normal probability density
y <- dnorm(z,6,1)
# scale it to sum to one
y <- y/sum(y)

# fitness function
w <- function(int,slp,z) int+slp*z 

# write a function to map x to x': Gaussian
gauss.fun <- function(x,y,params.m,params.v){
  mux <- params.m[1]+params.m[2]*x
  sigmax2 <- params.v[1]
  sigmax <- sqrt(sigmax2)
  fac1 <- sqrt(2*pi)*sigmax
  fac2 <- ((y-mux)^2)/(2*sigmax2)
  return(exp(-fac2)/fac1)
}

# num bins
n <- 500

# num iterations
n.gens <- 100

# midpoints
z <- seq(0,30,length.out=n)

# results vector
res.p <- res <- array(NA,c(n,n.gens))

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

# starting condition
res[,1] <- dnorm(z,8,1)
res[,1] <- res[,1]/sum(res[,1])
res.p[,1] <- res[,1]
# Calculate fitness
fit <- w(0.2,0.1,z)
# run the simulation
for (i in 2:n.gens){
  # moments of distribution post-selection
  temp <- moments.fun(fit*res[,i-1],z)
  # moments of distribution pre-selection
  temp2 <- moments.fun(res[,i-1],z)
  # iterate forward the clonal model
  res.p[,i] <- fit*res.p[,i-1]
  # moments of the clonal model
  temp.sel <- moments.fun(res.p[,i],z)
  # Gaussian transition approximation
  H <- t(outer(z,z,gauss.fun,c(temp[1]/2,0.5),0.75*temp[2]))
  H <- H/matrix(as.vector(apply(H,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # iterate forward the Gaussian approximation model
  res[,i] <- H%*%(fit*res[,i-1])
}

# now do the convolution
m <- 20000 # the larger the size of the matrix the more accurate the sampling of the convolution back to the original length becomes
y <- seq(0,30,length.out=m)
# length following convolution 1
y1 <- seq(0,30,length.out=2*m-1)
# length following convolution 2
y2 <- seq(0,30,length.out=4*m-3)
res.c <- array(NA,c(m,n.gens))
res.c[,1] <- dnorm(y,8,1)
res.c[,1] <- res.c[,1]/sum(res.c[,1])
fit <- w(0.2,0.1,y)

for (i in 2:n.gens){
  # parental distribution post selection (males and females have same demography)
  res.c[,i] <- fit*res.c[,i-1]
  size.of.dist <- sum(res.c[,i])
  scaled.dist <- res.c[,i]/sum(res.c[,i])
  # do first convolution - the result is 2n-1 long
  conv.1 <- convolve(scaled.dist,rev(scaled.dist),type='open')  
  moments.of.dist <- moments.fun(scaled.dist,y)
  mean.position <- findInterval(moments.of.dist[1],y1)
  # now construct the second distribution
  norm <- dnorm(y1,moments.of.dist[1],sqrt(moments.of.dist[2]/2))
  norm <- norm/sum(norm)
  conv.2 <- convolve(conv.1,rev(norm),type='open')
  ii <- seq.int(mean.position,2*m-1+mean.position,by=2)
  res.c[,i] <- conv.2[ii]/sum(conv.2[ii])
  res.c[,i] <- res.c[,i]*size.of.dist
}  
# calculate the population grwoth rates
Ns <- apply(res,2,sum)
R0 <- Ns[-1]/Ns[-length(Ns)]
Ns.sel <- apply(res.p,2,sum)
R0.sel <- Ns.sel[-1]/Ns.sel[-length(Ns.sel)]
Ns.selc <- apply(res.c,2,sum)
R0.selc <- Ns.selc[-1]/Ns.selc[-length(Ns.selc)]
# plot the population growth rates for the clonal, sexual and approximate models 
par(mfg=c(1,1))
plot(1:(n.gens-1),R0,type='l',ylim=c(min(R0.sel),max(R0.sel)),xlab='',ylab='',col='red')
lines(1:(n.gens-1),R0.sel,col='black')
lines(1:(n.gens-1),R0.selc,col='blue')
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'R0')
text(10,1.58,'(a)',cex=1.5)
leg.txt <- c('Clonal','Convolution','Gaussian')
leg.col <- c('black','red','blue')
legend('bottomright',col=leg.col,legend=leg.txt,lty=c(1,1,1),bty='n')

# now calculate the moments of the three distributions
kurtosis.cp <- skews.cp <- means.cp <- vars.cp <- kurtosis.c <- skews.c <- means.c <- vars.c <- kurtosis <- kurtosis.s <- skews <-skews.s <- means <- vars <- means.s <- vars.s <- rep(NA,n.gens)
for (i in 1:n.gens){
  temp <- moments.fun(res[,i],z)
  means[i] <- temp[1]
  vars[i] <- temp[2]
  skews[i] <- temp[3]
  kurtosis[i] <- temp[4]
  temp <- moments.fun(res.p[,i],z)
  means.s[i] <- temp[1]
  vars.s[i] <- temp[2]
  skews.s[i] <- temp[3]
  kurtosis.s[i] <- temp[4]
  temp <- moments.fun(res.c[,i],y)
  means.c[i] <- temp[1]
  vars.c[i] <- temp[2]
  skews.c[i] <- temp[3]
  kurtosis.c[i] <- temp[4]
}
# plot the various moments and the dynamics of the distribution
par(mfg=c(1,2))
plot(1:n.gens,means.s,type='l',xlab='',ylab='',col='black')
lines(1:n.gens,means,col='red')
lines(1:n.gens,means.c,col='blue')
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'Means')
text(10,13.8,'(b)',cex=1.5)
temp <- c(vars,vars.s)
par(mfg=c(1,3))
plot(1:n.gens,vars.s,type='l',ylim=c(min(temp),max(temp)),xlab='',ylab='')
lines(1:n.gens,vars,col='red')
lines(1:n.gens,vars.c,col='blue')
mtext(side=1,line=2.5,'Time')
mtext(side=2,line=2.5,'Variance')
text(90,0.99,'(c)',cex=1.5)

par(mfg=c(2,1))
plot(z,res[,1]/sum(res[,1]),type='l',ylim=c(0,0.03),xlab='',ylab='',lwd=2,col='green',xlim=c(0,35))
lines(z,res[,n.gens]/sum(res[,n.gens]),col='blue',lwd=2)
lines(z,res.p[,n.gens]/sum(res.p[,n.gens]),col='black',lwd=2)
mtext(side=1,line=2.5,'Breeding value')
mtext(side=2,line=2.5,'Density')
text(2,0.029,'(d)',cex=1.5)
leg.txt <- c('Initial','Clonal','Gaussian')
leg.col <- c('green','black','blue')
leg.lty <- c(1,1,1)
legend('topright',legend=leg.txt,col=leg.col,lty=leg.lty,bty='n')

# finally give an example of how dominance would influence an inheritance function
x.vals <- 0:12
y.vals.im <- 5+0.5*x.vals
y.vals.dv <- y.vals.im-0.5
par(mfg=c(3,1))
plot(x.vals,y.vals.im,type='l',xlab='',ylab='',xlim=c(0,12),ylim=c(0,12))
lines(x.vals,y.vals.dv,col='red')
mtext(side=1,line=2.5,'Parental genotype')
mtext(side=2,line=2.5,'Offspring genotype')
leg.txt <- c('No dominance var','Dominance var')
leg.col <- c('black','red')
leg.lty <- c(1,1)
legend('topleft',legend=leg.txt,col=leg.col,lty=leg.lty,bty='n')
text(11.4,1,'(g)',cex=1.5)

x.vals <- seq(0,12,length.out=100)
y.vals <- dnorm(x.vals,5,1)*10
lines(x.vals,y.vals,col='blue',lwd=2)

x <- c(5,5)
y <- c(0,7.5)
lines(x,y,lty=3,lwd=2)

x <- c(0,5)
y <- c(7.5,7.5)
lines(x,y,lty=3,lwd=2)

x <- c(0,5)
y <- c(7,7)
lines(x,y,col='red',lty=3,lwd=2)


