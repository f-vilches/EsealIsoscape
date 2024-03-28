# Script to set up xtracto, written by Dana Briscoe

###RUN RSTUDIO AS ADMINISTRATOR

rm(list = ls(all = TRUE)) #clear workspace

#read in R packages 
library("ggplot2")
library("ggfortify")	
library("lubridate")
library("mapdata")
library("RColorBrewer")
library("reshape2")
library("xts")
library("xtractomatic")
require("httr")
require("ncdf4") 
require("sp")

#set directories and read in tag data
path=("~/Grants and Fellowships/2016 NSF GRIP/")		
data.in.dir=paste(path, "Xtracto Code/", sep="")
data.out.dir=paste(path, "Xtracto Code/", sep="")

setwd(data.in.dir) # set working directory
	
tag = read.csv("SampleData_2008042.csv") # load csv file
head(tag) # check out data frame

# rename variables for consistency
tag$lat = tag$Latitude
tag$lon = tag$Longitude

# give this data frame a ptt/tag id
tag$ptt = rep('2007046', length(tag$Date))

# convert MATLAB date to R date using the following function
#'http://lukemiller.org/index.php/2011/02/converting-matlab-and-r-date-and-time-values/'
matlab2POS = function(x, timez = "UTC") {
  days = x - 719529 	# 719529 = days from 1-1-0000 to 1-1-1970
  secs = days * 86400 # 86400 seconds in a day
  
  return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
  		   tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
  	     tz = 'UTC', usetz = FALSE), tz = timez))
}
		 
tag$dt = matlab2POS(tag$Date) # format tag dates using a loop
head(tag)	# check results

# 'Pare down (reduce)' the dataset to 1-hit per day (first occurrence of date/time)

# remove the H-M-S component of the date-time. Store values in new column: 'tag$date'
tag$dtime = as.POSIXct(strptime(as.character(tag$dt), "%Y-%m-%d ", tz="America/Los_Angeles"))
head(tag$dtime)

# remove duplicate dates
dupes = which(duplicated(tag$dtime))		
if (length(dupes)!=0){
  tag = tag[-dupes,]
}

tag = subset(tag, select = c("ptt","lat", "lon", "dt",'FPT' ))# get rid of duplicate col/col names

## SET UP XTRACTO ##################################

datasetnames = 'ncdcOwDlyStrs' # list datasetnames    
#ncdcOisst2Agg - sst optimum interpolation, avhrr: 1-day, 0.25 spatial res - WORKING
  #mean for water masses, SD for fronts
#erdTAssh1day - Sea Surface height - WORKING
    #SD for eddies
#mhchla8day - WORKING, but deprecated
#etopo360 - TOPOGRAPHY
#ncdcOwDlyStrs - windstress
    #upwelling

###PRODUCTS TO FIGURE OUT
#erdMWchla8day # LONGITUDE RANGE NOT SUFFICIENT
#erdMBchla8day #lONGITUDE SUFFICIENT, but Dana said not to use?
    #overall productivity


ptts=unique(tag$ptt) # get unique ptt ids (if multiple)

# start xtracto loop
for (k in 1:length(ptts)){

  xtr=tag[tag$ptt==ptts[k],]
  xtr$id = xtr$ptt  

  # Assign variables
  ptt = unique(xtr$id)  # get unique tag ids
  print(ptt)
  tpos1=xtr$dt
  tpos=as.character(tpos1); head(tpos)
  xpos=xtr$lon
  ypos=xtr$lat

 	# ENVvars Loop -------------------------------
      
	 allvars=rep(NA,length(xpos))

   for(d in 1:length(datasetnames)){
       
        # xtract!	
      	dtype = datasetnames[d]
      	print(paste('now running env: ', dtype, sep=""))
        extract_temp =xtracto(xpos,ypos,tpos,dtype, xlen=.1,ylen=.1)

    	  extract_cols=cbind(extract_temp$mean,extract_temp$stdev,extract_temp$median)
 		    colnames(extract_cols)=c(paste(datasetnames[d],"_mean",sep=""),paste(datasetnames[d],"_stdev",sep=""),paste(datasetnames[d],"_med",sep=""))
 
        # bind extracted env data to original track data frame
        dframet=cbind(xtr,extract_cols)	

        # bind extracted env data to previous extracted data... for each iteration of d
        allvars=cbind(allvars,extract_cols); head(allvars)

      }

    allvars=as.data.frame(allvars[,-1])	# removes dummy column
    dframe=cbind(xtr,allvars)				# binds all extracted data to original track data frame

## Save individual track/xtracto file 
setwd(data.in.dir)
out.file = paste('xtracto_track_',ptt,'.csv', sep="")
print(paste('now saving track: ', ptt, sep=""))

write.csv(dframe, out.file)

} 


#READ IN RESULTS & MAKE BOXPLOTS #####################################

dat=read.csv('xtracto_track_2007046.csv')

#NEED TO ADD FILTER HERE
SST_mean=dat$ncdcOisst2Agg_mean
SST_sd=dat$ncdcOisst2Agg_stdev
Chla=dat$mhchla8day_mean
SSH_sd=dat$erdTAssh1day_stdev

Treatment=as.factor(dat$FPT) #0 is not FPT, 1 is FPT

boxplot(SST_mean~Treatment);fit=aov(SST_mean~Treatment);summary(fit)
boxplot(SST_sd~Treatment);fit=aov(SST_sd~Treatment);summary(fit)
boxplot(SSH_sd~Treatment);fit=aov(SSH_sd~Treatment);summary(fit)
boxplot(Chla~Treatment);fit=aov(Chla~Treatment);summary(fit)
