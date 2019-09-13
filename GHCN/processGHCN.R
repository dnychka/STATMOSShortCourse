library(fields)

# parameters for the analysis
startYear <- 1980
endYear   <- 2016

minYearsData <- 5

# Data downloaded from NOAA, repositiory avaiable at:
# https://www1.ncdc.noaa.gov/pub/data/ghcn/v4/
# The qcf data is quality controlled and homogenized
raw <- readLines('Data/ghcnm.tavg.v4.0.1.20190826.qcf.dat')

stationID   <- rep(NA, length(raw))
stationData <- matrix(NA,nrow=length(raw),ncol=13)
colnames(stationData) <- c('year', 1:12)

# function to read a line
parseDataString <- function(dataString){
	sapply(seq(from=1, to=nchar(dataString), by=8), 
		function(i) as.numeric(substr(dataString, i, i+4)))
}

# loop through the dataset
for(i in 1:length(raw)){
	temp <- raw[i]
	stationID[i] <- substr(temp,1,11)
	year    <- as.numeric(substr(temp,12,15))
	data <- parseDataString(substr(temp,20,nchar(temp)))

	stationData[i,] <- c(year,data)
}


# SUBSET in time given start year and end year

dataYears <- (startYear):endYear
nYears <- length(dataYears)

goodYearInds <- which(stationData[,1] %in% dataYears)
stationDataTemporal <- stationData[goodYearInds,]
stationIDTemporal   <- stationID[goodYearInds]

# Load in the station info
infoRaw <- readLines('Data/ghcnm.tavg.v4.0.1.20190826.qcf.inv')

infoID <- rep(NA, length(infoRaw))
locMat <- matrix(NA, nrow=length(infoRaw), ncol=2)
elev   <- rep(NA, length(infoRaw))

for(i in 1:length(infoRaw)){
	temp <- infoRaw[i]
	infoID[i] <- substr(temp, 1, 11)

	lat <- as.numeric(substr(temp, 13, 20))
	lon <- as.numeric(substr(temp, 22, 30))

	elev[i] <- as.numeric(substr(temp, 32, 37))
	locMat[i,] <- c(lon,lat)
}

# SUBSET the info in space, resricting to only the US stations

fipsCode <- substr(infoID,1,2)

# pull only the us values
usInfoInds <- which(fipsCode == 'US')

# crude removal of AK, HI
finalInds <- usInfoInds[locMat[usInfoInds,1] > -130]

conusID <- infoID[finalInds]
conusLoc <- locMat[finalInds,]
conusElev <- elev[finalInds]

infoFrame <- data.frame(id=conusID,lon=conusLoc[,1],lat=conusLoc[,2],elev=conusElev)


# SUBSET the station records given the station ids
# use to build the record matrix where row is station and column is time

noDataRows <- c()

tempRecordMat <- matrix(NA, nrow=length(conusID),ncol=12*nYears)


# crude loop to pull out CONUS data
for(i in 1:length(conusID)){
	tempID <- conusID[i]


	dataInds <- which(stationIDTemporal == tempID)

	subData <- stationDataTemporal[dataInds,]

	if(is.null(nrow(subData))){
		subData <- matrix(subData,nrow=1,ncol=13)
	}

	if(nrow(subData) < 1){
		noDataRows <- c(noDataRows,i)
	}

	tempRecord <- rep(NA, nYears)
	for(j in 1:nYears){

	
		yearInds <- which(subData[,1] == dataYears[j])

		if(length(yearInds) < 1){
			tempYear <- rep(NA, 12)
		} else{
			tempYear <- subData[yearInds,2:13]
		}

		tempYearClean <- tempYear/100

		# crude removal of missing values
		tempYearClean[tempYearClean < -50] <- NA


		tempRecord[(1 + 12*(j-1)):(12*j)] <- tempYearClean
	}

	tempRecordMat[i,] <- tempRecord

}



# remove the rows from the ID frame and the record matrix
infoFrameClean <- infoFrame[-noDataRows,]
tempMatClean <- tempRecordMat[-noDataRows,]

# Do one final screen for stations that have less than 5 years of data
naCount <- apply(tempMatClean,1,function(x) sum(is.na(x)))

cutoff <- 12*(nYears - minYearsData)

sparseYears <- which(naCount > cutoff)

# deal with crazy elevations
badElev <- which(infoFrameClean[,4] > 4000)
infoFrameClean[badElev,4] <- NA

# build the final objects
infoFrameFinal <- infoFrameClean[-sparseYears,]
tempMatFinal   <- tempMatClean[-sparseYears,]

# rename and save
stationInfo <- infoFrameFinal
tempDataMat <- tempMatFinal

save(stationInfo,tempDataMat,file='ghcnConusClean_1980_2016.Rda')

# visualize a year to make sure we ok

zr <- range(tempMatFinal[,1:12],na.rm=TRUE)

pdf('testViz.pdf',18,12)
par(mfrow=c(3,4))
for(i in 1:12){
	jan80z <- tempMatFinal[,i]
	jan80xy <- infoFrameFinal[,2:3]

	badInds <- which(is.na(jan80z))


	pal <- tim.colors(256)[as.numeric(cut(c(jan80z[-badInds],zr),breaks = 256))]

	US()
	title(main=paste(month.name[i],'1980'))
	points(jan80xy[-badInds,],col=pal)
}
dev.off()


# vizualize elevation

pdf('elevation.pdf',10,6)

jan80z <- infoFrameFinal[,4]
jan80xy <- infoFrameFinal[,2:3]

badInds <- which(is.na(jan80z))


pal <- tim.colors(256)[as.numeric(cut(jan80z[-badInds],breaks = 256))]

US()
title(main="Elevation")
points(jan80xy[-badInds,],col=pal)
dev.off()
