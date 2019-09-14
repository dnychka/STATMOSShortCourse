################################################################################################
## How convolutions (and kernel smoothing) using fft work in R
################################################################################################

##
## 1D
##

n <- 20
dat <- 1:n

## Example indicator kernels
kern <- rep(0,n)
kern[3] <- 1
kern
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n)
cbind(dat,conv)

kern <- rep(0,n)
kern[5] <- 1
kern
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n)
cbind(dat,conv)

## Example smoothing kernels: fill in the first (2*m+1) entries
m <- 1
kern <- rep(0,n)
kern[1:(2*m+1)] <- 1/(2*m+1) # for smoothing we want the kernel to sum to 1
kern
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n)
cbind(dat,conv)

m <- 4
kern <- rep(0,n)
kern[1:(2*m+1)] <- 1/(2*m+1) # for smoothing we want the kernel to sum to 1
kern
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n)
cbind(dat,conv)
mean((8-m):(8+m))
# note: the smoother of the Kth point is indexed at (K+m)
# note: the kernel "wraps" around the data vector

## Example smoothing kernels: fill in middle values with 1/(2*m+1)
center <- n-floor(n/2)
center
m <- 1
kern <- rep(0,n)
kern[(center-m):(center+m)] <- 1/(2*m+1) # for smoothing we want the kernel to sum to 1
cbind(kern)
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n)
cbind(dat,conv)

m <- 3
kern <- rep(0,n)
kern[(center-m):(center+m)] <- 1/(2*m+1) # for smoothing we want the kernel to sum to 1
cbind(kern)
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n)
cbind(dat,conv)
mean(c(18,19,20,1,2,3,4))
# note: the smoother of the Kth point is indexed at K+(center-1)

rm(list=ls())

##
## 2D
##

n <- 6
dat <- matrix(1:(n^2),n,n)

## Example indicator kernels
kern <- matrix(0,n,n)
kern[4,4] <- 1
kern
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n^2) # note n^2 here since n is marginal dimension
dat
conv

## Example smoothing kernels centerd at (center,center) of half-width m in each axial direction
center <- 4
m <- 1
kern <- matrix(0,n,n)
kern[(center-m):(center+m),(center-m):(center+m)] <- 1/(2*m+1)^2 # sums to one
kern
conv <- Re(fft(fft(kern)*fft(dat),inverse=TRUE)/n^2) # note n^2 here since n is marginal dimension
dat
conv
mean(c(1,2,7,8,31,32,6,12,36)) # (1,1) entry
mean(c(1:3,7:9,13:15)) # (2,2) entry
# note: the smoother of the (i,j)th point is indexed at (i+center-1,j+center-1)
# or equivalently, the (1,1) entry's smoother is at (center,center)

## Below, the (0,0) frequency in the DFT will be at the (1,1) entry, so when
## we do a kernel convolution this will move to the (center,center) position,
## which can make plots look nice, but it's also a pain to keep track of

################################################################################################
## 1D coherence example using only cos terms
################################################################################################

library(fields)
library(MASS)

##
## Set up basis functions
##

n <- 500 # number of equally-spaced spatial samples
nf <- n # number of frequencies

grd <- 1:n # observation grid
freqs <- 2*pi*c(0:(nf-1))/n

Phi <- matrix(nc=nf,nr=n)
for(l in 1:nf){
  Phi[,l] <- cos(grd*freqs[l])
}

matplot(Phi[,1:10],type="l")

##
## Constant coherence at 0.9
##

## These will hold weighted cosines for process 1 and 2
z1 <- z2 <- matrix(nc=nf,nr=n)

set.seed(101)
for(l in 1:nf){
  tmp <- mvrnorm(n=1,mu=c(0,0),Sigma=matrix(c(1,0.9,0.9,1),2,2)/(0.1+freqs[l]^2)^3)
  z1[,l] <- Phi[,l] * tmp[1]
  z2[,l] <- Phi[,l] * tmp[2]
}

## Full simulation
Z1 <- apply(z1,1,sum)
Z2 <- apply(z2,1,sum)

## Weighted components
par(mfrow=c(2,4))
for(i in 1:8){
  plot(z1[,i],type="l",ylim=range(z1[,i],z2[,i]))
  lines(z2[,i],type="l",col="red")
}
dev.off()
plot(Z1,type="l")
lines(Z2,type="l",col="red")

## Periodograms and cross-periodogram
p1 <- Mod(fft(Z1))^2/n
p2 <- Mod(fft(Z2))^2/n
p12 <- fft(Z1)*fft(Z2,inverse=TRUE)/n

## R calculates the DFT (via fft) for a vector of length n at frequencies 2*pi*(0:(n-1))/n
par(mfrow=c(2,1))
plot(freqs,p1,type="l",main="Periodogram for variable 1")
plot(freqs,p2,type="l",main="Periodogram for variable 2")
# note the symmetry: why is that?
# check Phi[,2] and Phi[,n]
# note: I am being lazy and not including the sin part

## Kernel smoothing two ways:
kern <- rep(0,n)
center <- n-floor(n/2)

m <- 15
kern[(center-m):(center+m)] <- 1/(2*m+1)
plot(kern,type="l")

kp1 <- Re(fft(fft(kern)*fft(p1),inverse=TRUE)/n)
kp2 <- Re(fft(fft(kern)*fft(p2),inverse=TRUE)/n)
kp12 <- Re(fft(fft(kern)*fft(p12),inverse=TRUE)/n)

coh <- sqrt(Mod(kp12)^2 / (kp1 * kp2))
# This has the effect of reindexing 
plot(freqs[1:(center+1)],coh[center:n],type="l",ylim=c(0,1),xlab="Frequency",ylab="Coherence")
abline(h=0.9)

##
## Decaying coherence like 0.9*exp(-5*frequency)
##

## These will hold weighted cosines for process 1 and 2
z1 <- z2 <- matrix(nc=nf,nr=n)

set.seed(102)
for(l in 1:nf){
  Sigma <- matrix(c(1,0.9*exp(-5*freqs[l]),0.9*exp(-5*freqs[l]),1),2,2) # coherence part
  Sigma <- Sigma / (0.1+freqs[l]^2)^3 # spectral density (variance) part
  tmp <- mvrnorm(n=1,mu=c(0,0),Sigma=Sigma)
  z1[,l] <- Phi[,l] * tmp[1]
  z2[,l] <- Phi[,l] * tmp[2]
}

## Full simulation
Z1 <- apply(z1,1,sum)
Z2 <- apply(z2,1,sum)

## Weighted components
par(mfrow=c(2,4))
for(i in 1:8){
  plot(z1[,i],type="l",ylim=range(z1[,i],z2[,i]))
  lines(z2[,i],type="l",col="red")
}
dev.off()
plot(Z1,type="l")
lines(Z2,type="l",col="red")

## Periodograms and cross-periodogram
p1 <- Mod(fft(Z1))^2/n
p2 <- Mod(fft(Z2))^2/n
p12 <- fft(Z1)*fft(Z2,inverse=TRUE)/n

## R calculates the DFT (via fft) for a vector of length n at frequencies 2*pi*(0:(n-1))/n
par(mfrow=c(2,1))
plot(freqs,p1,type="l",main="Periodogram for variable 1")
plot(freqs,p2,type="l",main="Periodogram for variable 2")
# note the symmetry

## Kernel smoothing two ways:
kern <- rep(0,n)
center <- n-floor(n/2)

m <- 10
kern[(center-m):(center+m)] <- 1/(2*m+1)
plot(kern,type="l")

kp1 <- Re(fft(fft(kern)*fft(p1),inverse=TRUE)/n)
kp2 <- Re(fft(fft(kern)*fft(p2),inverse=TRUE)/n)
kp12 <- Re(fft(fft(kern)*fft(p12),inverse=TRUE)/n)

coh <- sqrt(Mod(kp12)^2 / (kp1 * kp2))
# This has the effect of reindexing 
plot(freqs[1:(center+1)],coh[center:n],type="l",ylim=c(0,1))
lines(freqs[1:(center+1)],0.9*exp(-5*freqs[1:(center+1)]))

################################################################################################
## Spectral coherence study for a bivariate Matern simulation
################################################################################################

library(RandomFields)
library(fields)

## Theoretical spectral densities for the Matern
## - note we need common range in this definition!
## - the spectrum is radially symmetric, w = frequency length = ||(w.x,w.y)||
## - w = frequency
## - a = inverse range
## - nu11/nu22 = smoothnesses for process 1 and 2
## - nu12 = cross-smoothness
## - rho = marginal cross-correlation coefficient
## - d = dimension
spec.matern <- function(w,a,nu11,nu22,nu12,rho,d){
  f12 <- rho * (gamma(nu12+d/2) / gamma(nu12)) * (a^(2*nu12) / pi^(d/2)) /
    (a^2 + w^2)^(nu12 + d/2)
  f11 <- (gamma(nu11+d/2) / gamma(nu11)) * (a^(2*nu11) / pi^(d/2)) /
    (a^2 + w^2)^(nu11 + d/2)
  f22 <- (gamma(nu22+d/2) / gamma(nu22)) * (a^(2*nu22) / pi^(d/2)) /
    (a^2 + w^2)^(nu22 + d/2)
  return(list=list(f11=f11,f22=f22,f12=f12,coh=sqrt(f12^2/(f11*f22))))
}

##
## Bivariate simulations using RandomFields (exploits circulant embedding)
##

## Bivariate Matern setup
nu <- 0.5 # nu in RMbiwm is nu11, nu12, nu22
range <- 2 # actual range, not inverse range
rho <- 0.3 # not all correlation coefficients will be valid, try 0.7, e.g.
model <- RMbiwm(nu=c(nu,nu+1,nu), s=c(range,range,range), c=c(1,rho,1))
# biwm stands for bivariate Whittle-Matern

ns <- 200 # number of spatial samples in one axial direction
2*ns^2 # number of random variables we are simulating
x.seq <- y.seq <- 1:ns

x.mat <- matrix(x.seq,ns,ns)
y.mat <- t(matrix(y.seq,ns,ns))

set.seed(10)
sim <- RFsimulate(model,x.seq,y.seq)
sim1 <- matrix(sim$variable1,ns,ns)
sim2 <- matrix(sim$variable2,ns,ns)

par(mfrow=c(1,2))
image.plot(x.mat,y.mat,sim1,xlab="",ylab="",main="Process 1")
image.plot(x.mat,y.mat,sim2,xlab="",ylab="",main="Process 2")

##
## DFTs and manual smoothing
##

# Note: the way R organizes its DFT via fft is on frequencies
# outer product of (2*pi*0/n, 2*pi*1/n, ..., 2*pi*(n-1)/n)
# so the [1,1] entry should equal sum(dat)^2/ns^2
# (note this is particular to having the domain be square and having "n" be the marginal dimension)

## Periodograms: estimates of f_11, f_22 and f_12
p1 <- Mod(fft(sim1))^2/ns^2
p2 <- Mod(fft(sim2))^2/ns^2
p12 <- fft(sim1)*fft(sim2,inverse=TRUE)/ns^2

## Also note we *have* to do smoothing if we want coherence:
range( Mod(p12)^2 / (p1 * p2) )

## Kernel smooth the resulting periodogram and cross-periodograms with a box of width 2*m+1
kern.mat <- matrix(0,ns,ns)
m <- 3 # sure why not: don't want to undersmooth or oversmooth
c.x <- c.y <- ns-floor(ns/2)
kern.mat[(c.x-m):(c.x+m),(c.y-m):(c.y+m)] <- 1/(2*m+1)^2
image.plot(kern.mat)

## Kernel-smoothed periodograms for f_11, f_22 and f_12
kp1 <- Re(fft(fft(kern.mat)*fft(p1),inverse=TRUE)/ns^2)
kp2 <- Re(fft(fft(kern.mat)*fft(p2),inverse=TRUE)/ns^2)
kp12 <- Re(fft(fft(kern.mat)*fft(p12),inverse=TRUE)/ns^2)

coh <- sqrt(Mod(kp12)^2 / (kp1 * kp2))
# This has the effect of reindexing the frequencies so that the (0,0) frequency is
# *not* at [1,1], e.g., but now it is at kp1[ns-floor(ns/2),ns-floor(ns/2)]
# sanity check: compare these values against sum(sim1)^2/ns^2
p1[1,1]
kp1[c.x,c.y] # IF m = 0 above so that no kernel smoothing is done
sum(sim1)^2/ns^2

##
## 2D plots and keeping track of frequencies
##

# First set up the frequency indexing
# R.freqs = R's DFT indexing
# kern.freqs = reindexed frequencies after kernel smoothing
R.freqs <- 2*pi*as.matrix(expand.grid((0:(ns-1))/ns,(0:(ns-1))/ns))
R.freqs.x <- matrix(R.freqs[,1],ns,ns)
R.freqs.y <- matrix(R.freqs[,2],ns,ns)
# The reindexed frequencies
kern.freqs <- 2*pi*as.matrix(expand.grid((-floor((ns-1)/2)):floor(ns/2),(-floor((ns-1)/2)):floor(ns/2)))/ns
kern.freqs.x <- matrix(kern.freqs[,1],ns,ns)
kern.freqs.y <- matrix(kern.freqs[,2],ns,ns)

# Now do the plotting
par(mfrow=c(2,3))
image.plot(R.freqs.x,R.freqs.y,p1)
image.plot(R.freqs.x,R.freqs.y,p2)
plot.new() # note we can't plot p12 easily since it's complex
image.plot(kern.freqs.x,kern.freqs.y,kp1)
image.plot(kern.freqs.x,kern.freqs.y,kp2)
image.plot(kern.freqs.x,kern.freqs.y,coh,zlim=c(0,1))

##
## Check against theoretical squared coherences
##

## The periodograms and estimated coherence are over two dimenions
## But, they are radially symmetric, so just point in one angle and walk along radius (frequency length)
freqs <- (0:floor(ns/2))/ns
sm <- spec.matern(w=2*pi*freqs,a=1/range,nu11=nu,nu12=nu+1,nu22=nu,rho=rho,d=2)

par(mfrow=c(2,2))
plot(c(kp1[c.x,(c.y:ns)]/(2*pi)^2)~freqs,type="l",ylim=c(0,max(sm$f11)),xlab="Frequency",ylab="P_11",
  main="Density for variable 1")
lines(sm$f11~freqs,col="red")
plot(c(kp2[c.x,(c.y:ns)]/(2*pi)^2)~freqs,type="l",ylim=c(0,max(sm$f11)),xlab="Frequency",ylab="P_22",
  main="Density for variable 2")
lines(sm$f22~freqs,col="red")
plot(c(kp12[c.x,(c.y:ns)]/(2*pi)^2)~freqs,type="l",ylim=c(0,max(sm$f11)),xlab="Frequency",ylab="P_12",
  main="Cross-spectral density (modulus of)")
lines(sm$f12~freqs,col="red")
plot(c(coh[c.x,(c.y:ns)])~freqs,type="l",ylim=c(0,1),xlab="Frequency",ylab="Coherence",main="Coherence")
lines(sm$coh~freqs,col="red")

################################################################################################
## Just for fun...
## 2D stationary simulation: Rigged with band of coherence
################################################################################################

library(fields)

##
## Setup
##

## Spatial grid
ns <- 100 # per axis
grd <- as.matrix(expand.grid(1:ns,1:ns))

## Frequency grid
freqs <- 2*pi*(-floor((ns-1)/2)):floor(ns/2)/ns
omegas <- as.matrix(expand.grid(freqs,freqs))
omegas.length <- sqrt(rowSums(omegas^2))

sw <- grd %*% t(omegas) # (i,j)th entry is s_i^T omega_j

## Coherence matrix for frequencies that will be correlated
Sigma.cor <- matrix(c(1,0.9,0.9,1),2,2)
Sigma.uncor <- diag(1,2)

## Variance(amplitude(omega_i)) = 1 / (inverse-range + ||omega_i||^2)^power
powers <- c(1,1)
sds.1 <- sqrt(1 / (0.01 + omegas.length^2)^powers[1])
sds.2 <- sqrt(1 / (0.01 + omegas.length^2)^powers[2])
sds <- cbind(sds.1,sds.2)

## Find indices for which we'll correlate coefficients
these <- which(abs(omegas.length-1) < 0.5) # ||omega_i|| approximately one +- 0.5

## Generate random amplitudes
set.seed(58)
amplitudes <- array(dim=c(ns^2,2))
for(i in 1:ns^2){
  if(i %in% these){
    amplitudes[i,] <- t(chol(sds[i,] %*% t(sds[i,]) * Sigma.cor)) %*% rnorm(2)
  }else{
    amplitudes[i,] <- t(chol(sds[i,] %*% t(sds[i,]) * Sigma.uncor)) %*% rnorm(2)
  }
}
# note that doing t(chol(...)) each time is the wrong way, but it is instructive

## Generate the bivariate simulation
csw <- cos(sw)
sim1 <- matrix(colSums(amplitudes[,1]*t(csw)),ns,ns)
sim2 <- matrix(colSums(amplitudes[,2]*t(csw)),ns,ns)

# plot the random fields -- is this a reasonable model?
par(mfrow=c(1,2))
quilt.plot(grd,sim1,nx=ns,ny=ns)
quilt.plot(grd,sim2,nx=ns,ny=ns)

## Periodograms and kernel-smoothed periodograms
p1 <- Mod(fft(sim1))^2/ns^2
p2 <- Mod(fft(sim2))^2/ns^2
p12 <- fft(sim1)*fft(sim2,inverse=TRUE)/ns^2

kern.mat <- matrix(0,ns,ns)
m <- 2
c.x <- c.y <- ns-floor(ns/2)
kern.mat[(c.x-m):(c.x+m),(c.y-m):(c.y+m)] <- 1/(2*m+1)^2

kp1 <- Re(fft(fft(kern.mat)*fft(p1),inverse=TRUE)/ns^2)
kp2 <- Re(fft(fft(kern.mat)*fft(p2),inverse=TRUE)/ns^2)
kp12 <- Re(fft(fft(kern.mat)*fft(p12),inverse=TRUE)/ns^2)

coh <- sqrt(Mod(kp12)^2 / (kp1 * kp2))

par(mfrow=c(1,3))
quilt.plot(omegas,kp1,xlab="x frequency",ylab="y frequency",main="Periodogram for variable 1")
quilt.plot(omegas,kp2,xlab="x frequency",ylab="y frequency",main="Periodogram for variable 2")
quilt.plot(omegas,coh,zlim=c(0,1),xlab="x frequency",ylab="y frequency",main="Coherence")
# identifies the frequencies for which the coefficients were correlated

# plot as function of frequency length
# recall correlation was only at frequencies with length in (0.5,1.5)
plot(c(coh)~omegas.length,ylim=c(0,1),xlab="Frequency length",ylab="Coherence")
abline(v=1,col="grey")
abline(v=0.5,col="red")
abline(v=1.5,col="red")

# note that marginal correlation is pretty uninformative
cor(c(sim1),c(sim2))

##
## Note that coherence tells us way more information (in this case!) than cross-covariance
##

## ONLY execute if ns small enough (e.g., 64)
emp.cross.cor <- crossCoVGram(loc1=grd,loc2=grd,y1=c(sim1),y2=c(sim2),lon.lat=FALSE,dmax=50,
  type="cross-correlogram")
plot(emp.cross.cor)
abline(h=0)
plot(emp.cross.cor,N=50)
abline(h=0)



