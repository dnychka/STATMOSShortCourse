
#########################################
### Lattice vs spatialProcess demo
#########################################

library(LatticeKrig)
data(CO2)

s<- CO2$lon.lat
z<- CO2$y
dim( s)
bubblePlot( s,z, highlight=FALSE, size=.4)
world( add=TRUE, col="magenta")


# 1884 locations
ind2<- which(
       s[,1]>= -120 & s[,1] <= -60 &
       s[,2]>=  -10 & s[,2] <=  55
        )

system.time(
  fit2<- spatialProcess(s[ind2,], z[ind2], cov.function="Tps.cov" )
)
# for n=1884 get about 7 seconds
# a conservative  lower bound   for the full data set is 
# 3 hours ...

surface( fit2)
world( add=TRUE, col="magenta")
title("subset of CO2 data")


# approx thin plate spline fit using fixed rank Kriging
system.time(
  fit4<- LatticeKrig(s, z, a.wght = 4.01 )
)
# for ~27K locations get about  35 seconds

# summary of fit
fit4

surface( fit4)
world( add=TRUE, col="magenta")


# fit model using cylindrical geometry.



alpha<- fit4$LKinfo$alpha
LKInfo<- LKrigSetup( s, startingLevel=4 , nlevel=3,
                                              a.wght=1.05,
                          alpha=alpha,
                                              LKGeometry="LKSphere" )
# list model specification
LKInfo

# system.time(
#   fit5<- LatticeKrig(s, z, LKinfo=LKinfo)
# )

# this takes abot 2 minutes and is an exact  spherical geometry. 


How different are the fits using the usual 2D model and the spherical one?
  
Do you see e

