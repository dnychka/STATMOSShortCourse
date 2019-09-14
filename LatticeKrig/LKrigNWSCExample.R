load("NWSC2.rda")
library( LatticeKrig)
# LatticeKrig and also spatialProcess do not handle 
# spatial replicates. Replace with average. 

out<- Krig.replicates(x=cbind(NWSC$Otemp, NWSC$RH ) ,
                      y=NWSC$Mpower)
x<- out$xM
y<- out$yM

# EDA plot
#pdf("pix/figMpower.pdf", width=7, height=6)
fields.style()
par(mar=c(5,4,3,1))
quilt.plot( cbind(NWSC$RH, NWSC$Otemp ) , NWSC$Mpower, nrow=100, ncol=100,
                  xlab="Outside Relative Humidity", ylab="Outside Temperature",
                  main= "Mechanical Systems Power Use,  October 2012",
                  legend.mar=8.1,
                  legend.args=list( text="kW", cex=1.2, side=4, 
                                    line=3))
#dev.off()

# standard computation using thin plate spline  takes about 
# 6 seconds
system.time(
  fit0 <- Tps(x,y)
)

# spatial process estimate also finding correlation range
system.time(
  fit0 <- spatialProcess(x,y)
)

# takes about 100 seconds, 10 seconds if a.wght fixed
system.time(
  fit<- LatticeKrig( x,y, NC=10, nlevel=4, a.wght=8.4 )
)
system.time(
fit<- LatticeKrig( x,y, NC=10, nlevel=4, findAwght=TRUE )
)

print( fit)
quilt.plot( x, fit$residuals)

surface( fit, xlab="Temp", ylab="RH")
# dev.copy2pdf( file="pix/NWSCfit.pdf", width=6, height=4)

simFit<- LKrig.sim.conditional( fit,  M=50)
image.plot( as.surface( simFit$x.grid, simFit$SE))


