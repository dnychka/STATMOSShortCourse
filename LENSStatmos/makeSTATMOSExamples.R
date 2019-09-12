setwd("/Users/nychka2/Data/LENS")
library( fields)
library( ncdf4)
source("subsetLENS.R")

# default region ~~ North America
# grid boxes I1 = 155, I2 = 250, J1 = 110, J2 = 160,
# Only extract 10 members

# surface Temperature 
LENSStatmosTS<- subsetLENS(variableName="TS",
                  nMember=10)
cat("writing to file", fill=TRUE)
save( LENSStatmosTS, 
      file="data/LENSStatmosTS.rda")

load("data/LENSStatmosTS.rda")
print( names( LENSStatmosTS) )

# plot of checking
pdf("data/STATMOSCheckTS.pdf",height=6, width=8)
set.panel(2,2)
image.plot(LENSStatmosTS$lon, LENSStatmosTS$lat, 
           LENSStatmosTS$field[,,7,10],
           zlim=c(273,315))
map("world2", add=TRUE, col="grey70")
title( "july 2006 10th member")
JJA<-  apply(LENSStatmosTS$field[,,6:8,],c(1,2), mean)
image.plot(LENSStatmosTS$lon, LENSStatmosTS$lat, JJA,
           zlim=c(273,315))
map("world2", add=TRUE, col="grey70")
title("Mean JJA")
timeTS<- LENSStatmosTS$year + (LENSStatmosTS$month-.5)/12
plot(timeTS[1:120], LENSStatmosTS$field[55,25,1:120,10],
     type="l", lwd=2, col="red")
title("grid box 55,25")
gridMeans<-  tapply( 
  LENSStatmosTS$field[55,25,,10],
  LENSStatmosTS$month,
  mean )
plot( 1:12, gridMeans, 
      type="l", lwd=2)
title("grid box 55,25")
dev.off()

# Large scale precip
LENSStatmosPRECL<- subsetLENS(variableName="PRECL",
                  nMember=10)
cat("writing to file", fill=TRUE)
save( LENSStatmosPRECL, 
      file="data/LENSStatmosPRECL.rda")
# Solar Insolation
LENSStatmosSOLIN<- subsetLENS(variableName="SOLIN",
                  nMember=10)
cat("writing to file", fill=TRUE)
save( LENSStatmosSOLIN, 
      file="data/LENSStatmosSOLIN.rda")
# 10 Meter Wind speed
LENSStatmosU10<- subsetLENS(variableName="U10",
                  nMember=10)
cat("writing to file", fill=TRUE)
save( LENSStatmosU10, 
      file="data/LENSStatmosU10.rda")

remove( list=ls())

### checks





