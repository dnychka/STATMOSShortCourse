library( fields)
data(ozone2)

s<- ozone2$lon.lat
y<- ozone2$y[16,]
# omit missing values to make the simulation below easier
ind<- !is.na(y)
s<- s[ind,]
y<- y[ind]

# only take out a constant in the fixed part of model
 obj<- spatialProcess( s,y, mKrig.args= list( m=1) )
 
 vObj<- vgram( s,y, N= 15)
 # add fitted variogram from spatialProcess
 
 dGrid<- seq( 0, 10, length.out=50)
 
 vFit<- obj$summary["tau"]^2 + 
   obj$summary["sigma2"]*( 1- Matern( dGrid/obj$summary["aRange"],
                                     smoothness=1.0 )
                        )
 plot(vObj)
 lines( dGrid, vFit, col="orange3", lwd=2)
 
 # now simulate data from the fitted model and examine variability in
 # the variograms
 vTest<- NULL
 MLETest<- NULL
 set.seed(222)
 ySim<- simSpatialData( obj, M=20)
 for( k in 1:20){
   cat( k, " ")
   
   vTemp<- vgram( s, ySim[,k], N= 15)
   vTest<- cbind( vTest, c(vTemp$stats[2,]) )
   
   objTest <- spatialProcess( s,ySim[,k], mKrig.args= list( m=1) )
   vFitTest<- objTest$summary["tau"]^2 + 
     objTest$summary["sigma2"]*( 1- Matern( dGrid/objTest$summary["aRange"],
                                        smoothness=1.0 )
     )
   MLETest<- cbind( MLETest, vFitTest)
   
 }
 library( scales)
 matplot(vObj$centers, vTest, type="b", lty=1, col="grey", pch=16 )
 matlines( dGrid, MLETest, col=alpha("magenta",.5), lty=1)
 lines( dGrid, vFit, col="black", lwd=2)
 
 
 matplot(vObj$centers, vTest, type="b", lty=1, col="grey", pch=16,
         xlim=c(0,4), ylim=c(0,2000))
 
 matlines( dGrid, MLETest, col=alpha("magenta",.5), lty=1)
 matlines( dGrid, vFit, col="black", lwd=2)
 
 