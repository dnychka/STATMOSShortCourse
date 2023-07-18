setwd("~/Dropbox/Home/Teaching/ISI Short Course June 1-3 Materials/Day 2 Materials")
library( LatticeKrig)
fields.style()
M<- 10
u<- seq( 0,1,length.out=M)
sGrid<- seq( 0,1,length.out= 100)

Omega<- exp( -rdist( u,u)/.2)
X<- Wendland2.2( rdist( sGrid, u)/.4)
matplot( sGrid, X, type="l")

covG<- X%*%Omega%*%t(X)
locs<- c( 10, 20, 50, 85)

library( viridisLite)


pdf("pix/FRKCovCoef.pdf", width=6, height =6)

image.plot( u, u , Omega, col=viridis(16), 
            cex=3)
#title("Covariance of coefficents",cex=10 )
dev.off()

pdf("pix/FRKPrecisionCoef.pdf", width=6, height =6)

image.plot( u, u , solve(Omega), col=viridis(16), 
            cex=3)
#title("Precision matrix  of coefficents",cex=10 )
dev.off()

pdf("pix/FRKCovImage.pdf", width=6, height =6)
fields.style()
image.plot( sGrid, sGrid, covG, col=viridis(256))
yline(sGrid[locs], col="grey", lwd=3)
#title("Covariance of g", cex=3)
dev.off()

pdf("pix/FRKCov.pdf", width=6, height =4)
fields.style()
par( mar=c(3,3,1,1))
matplot( sGrid, covG[,locs], type="l", lty=1, lwd=2,
         xlab= "s", ylab="covariance", col="orange3")

xline( sGrid[locs], col="grey", lwd=2)
dev.off()


pdf("pix/SAR1D.pdf", width=6, height =6)
fields.style()
B<- diag( 4, M)
B[ col(B) == row(B) +1]<- -1
B[ col(B) ==  row(B)-1]<- -1
temp<- B
temp[ temp==0]<- NA
image.plot(1:M, 1:M,  temp, col=viridis(256))
print( B)
dev.off()


pdf("pix/SAR1DPrecision.pdf", width=6, height =6)
fields.style()
B<- diag( 4, M)
B[ col(B) == row(B) +1]<- -1
B[ col(B) ==  row(B)-1]<- -1
Q<- t(B)%*%B
temp<- Q
temp[ temp==0]<- NA
image.plot(1:M, 1:M,  temp, col=viridis(256))
print( Q)
dev.off()

pdf("pix/SAR1DCovariance.pdf", width=6, height =6)
fields.style()
B<- diag( 4, M)
B[ col(B) == row(B) +1]<- -1
B[ col(B) ==  row(B)-1]<- -1
Q<- t(B)%*%B
image.plot(1:M, 1:M,  solve(Q), col=viridis(256))
print( solve(Q) )
dev.off()


set.seed(2232)
B<- chol( Omega)
c<- t(B)%*%rnorm( M)
g1<- X%*%c
fields.style()
plot( sGrid, g1, type="l")
c<- t(B)%*%rnorm( M)
g2<- X%*%c
c<- t(B)%*%rnorm( M)
g3<- X%*%c
matplot( sGrid, cbind( g1,g2,g3), type="l",
         lty=1, lwd=3)
#knot locations
xline( u, lwd=.5)
