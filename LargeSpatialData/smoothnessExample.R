library( fields)
data(ozone2)

s<- ozone2$lon.lat
y<- ozone2$y[16,]
# omit missing values to make the simulation below easier
ind<- !is.na(y)
s<- s[ind,]
y<- y[ind]

# only take out a constant in the fixed part of model
 obj <- spatialProcess( s,y, mKrig.args= list( m=1) )
 
 obj1<- spatialProcess( s,y, mKrig.args= list( m=1),
                        cov.function ="stationary.cov",
                       Covariance="Matern",
                       cov.params.start = list( 
                          smoothness=1.0, aRange=.5, lambda=.1 )
                      )
 
 obj2<- spatialProcess( s,y, mKrig.args= list( m=1),
                        parGrid = data.frame( smoothness = c(.5,.75,1.0,1.25) ),
                        cov.params.start = list( 
                        aRange = .5, lambda=.1), verbose=TRUE 
 )

 
 vFit<- obj$summary["tau"]^2 + 
   obj$summary["sigma2"]*( 1- Matern( dGrid/obj$summary["aRange"],
                                     smoothness=1.0 )
 
 ySim<- simSpatialData( obj, M=200)
 
 obj1<- spatialProcess( s,ySim, mKrig.args= list( m=1) )
 
 obj1B<- spatialProcess( s,ySim, mKrig.args= list( m=1),
   cov.params.start = list( aRange=.5, lambda=.1 )
   )
 

 obj1C<- spatialProcess( s,ySim, mKrig.args= list( m=1),
      cov.params.start = list(smoothness= 1.0, aRange=.5, lambda=.1 ),
      verbose=TRUE
 )
 
 
 
 obj2<- spatialProcess( s,ySim, mKrig.args= list( m=1),
                        cov.function ="stationary.cov",
                        Covariance="Matern",
                        cov.params.start = list( 
                           smoothness=1.0, aRange=.5, lambda=.1 )
 )
 
 