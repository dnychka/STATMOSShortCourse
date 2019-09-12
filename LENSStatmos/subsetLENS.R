
# TGCLDLWP,
# TS  (Radiative) surface temperature
# U10 10 meter Wind 
# PRECT	
# Total (convective and large-scale) precipitation rate (liq + ice)
# SOLIN Solar insolation

# get bounding box based on NARCCAP
#load( "/users/nychka2/Data/Projects/DEMworld/NARCCAP.elev.rda")
#lonRange<- range( NARCCAP.elev$x)
#lonRange360<- ifelse( lonRange <0,lonRange +360, lonRange) 
#latRange<- range( NARCCAP.elev$y)
#rm( "NARCCAP.elev")


subsetLENS<- function( variableName ='TS',
                       I1= 155,
                       I2 =  250,
                       J1 = 110,
                       J2 = 160,
                       nMember=NULL
                        ){
wildCardFile<-paste0("/Users/nychka2/Data/LENS/raw/*85*",
                         variableName,".200601-208012.nc")

fileNames<- system( paste0('ls ', wildCardFile ) , intern =TRUE)
N<- ifelse( is.null( nMember), length( fileNames),nMember)


# get common info from the first ensemble ncdf file
handle<- nc_open(fileNames[1]  )

latRaw <- ncvar_get( handle, 'lat')
lonRaw <- ncvar_get( handle, 'lon')

Longitude.0360<- lonRaw 
Longitude.trans <- ifelse(Longitude.0360<=180, Longitude.0360, Longitude.0360-360 ) 
index.sort <- order(Longitude.trans)
Longitude <- Longitude.trans[index.sort]
Latitude  <- latRaw

nx<- length( lonRaw)
ny<- length( latRaw)

#indices for NARCCAP bounding box
#I1<- which.max( lonRange360[1] <= lonRaw)
#I2<- which.min( lonRange360[2] >= lonRaw)

#J1<- which.max( latRange[1] <= latRaw)
#J2<- which.min( latRange[2] >= latRaw)

# these are smaller 
#I1<- 155
#I2<- 250
#J1<- 110
#J2<- 160

lonRegion<- Longitude.0360[I1:I2]
latRegion<- Latitude[J1:J2]

#handle<- nc_open(fileNames[member]  )
#var<- ncvar_get( handle,variableName,
#                  start=c(I1,            J1, 1 ),
#                  count=c(I2-I1 +1 ,J2-J1+1, 1 )
#                  )
#image.plot( lonRegion, latRegion, var)
#map( "world2", add=TRUE, col="magenta")



NTime <- length( date) 
date  <- ncvar_get( handle, 'date')
month<- substr( as.character( date), 5,6)
month<- as.numeric( month)
# assume that the month reported is actually the following month for the mean. 
# e.g. 2 is really January, 12 is November
month <- month - 1
# what was Jan is Dec.
month[ month == 0] <- 12
nTime<-  length( month)

if( nTime%%12!=0){
  stop("Not an even number of years")
}

year<- as.numeric(substr( as.character( date), 1,4))
# last raw year is jan of next year!
year<- c(year[1], year[1:(nTime-1)])
if( any(table(year)!=12) ) {
  stop("Not an even number of months")
}

out<- array( NA, c(I2-I1 +1 ,J2-J1+1, nTime,N ))
nX<- I2-I1 +1
nY<- J2-J1+1
cat("dim of field array", dim( out), fill=TRUE )
for( member in 1:N ){
  cat( ' ', fill=TRUE)
  cat( ' Working on  member ', member,  fill=TRUE)
  handle<- nc_open( fileNames[member]  )
  out[,,,member] <- ncvar_get( 
                  handle,
                  variableName,
                  start=c(I1, J1, 1 ),
                  count=c(nX , nY, nTime )
                  )
} 

out<- array( as.single( out), dim( out))

return( list( field= out,
              lon=lonRegion,
              lat=latRegion,
            month=month,
            year=year,
            nTime=nTime,
            nX= nX,
            nY=nY,
            nMember= N,
            J1=J1, J2=J2, I1=I1, I2=I2, 
            variableName= variableName,
            fileToken=wildCardFile,
            sourceCode = subsetLENS,
            creationDate=date()) 
        )

}
