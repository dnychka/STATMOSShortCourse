

fig1DApprox<- function(kappa= 25, A = TRUE){
B<- LKDiag( c( -1,2 + 1/kappa ,-1), 201, diags=c(-1,0,1))
Binv<- solve(B)
K0<- Binv%*%t(Binv)
Kd0<- sqrt(diag( K0))
K0 <-   t(K0/Kd0)/ Kd0

xg<- 100:-100
K<- stationary.cov( as.matrix( xg), 0,Covariance="Matern",
                    smoothness=1.0, theta = 6
                    )
K<- c( K)


if( A ){
  matplot(xg, cbind( K0[,101],K  ),
          xlab="distance", ylab="correlation", 
          type="l", lty=1 )
}
else{
  plot( K0[,101], K, pch=16, log="xy", col="grey")
  abline( 0,1, col="red",
          xlab="SAR Model", ylab="Matern")
}
#print( c( kappa, sqrt( mean((K0[,101]- K )^2)) )) 
}


