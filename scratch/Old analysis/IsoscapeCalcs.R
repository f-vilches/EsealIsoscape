######################################
### DRIFT PEAKS PLOT
######################################

setwd("~/Publications/In Prep/NESE Isoscape Project/Isoscape Raw Data")
dat=read.csv('IsoscapeDailyDriftRates.csv')

library(RColorBrewer)
uniqueseals=unique(dat$TOPPID)
Color=brewer.pal(length(uniqueseals),"BrBG")

png(file='DriftRatePeaks.png',bg="white",units="in", width=6, height=4, res=600)

for(i in 1:length(uniqueseals)){
  dat_i=dat[dat$TOPPID==uniqueseals[i],]
  if(i==1){
    plot(dat_i$DaysSinceJune1,dat_i$DriftRate,
         xlab="Days Since June 01",ylab="Daily Drift Rate (m/sec)",
         ylim=c(-0.02,0.1),type="l",col=Color[i],lwd=2)
  }else{
    par(new=TRUE)
    plot(dat_i$DaysSinceJune1,dat_i$DriftRate,
         xlab="",ylab="",xaxt="n",yaxt="n",
         ylim=c(-0.02,0.1),type="l",col=Color[i],lwd=2)
  }
}

abline(h=0)
dev.off()

##########################
### ISOSCAPE CALCULATIONS ALL POINTS
##########################
setwd("~/Publications/In Prep/NESE Isoscape Project/Isoscape Raw Data")
dat=read.csv('NESE Isoscape Isotope Data Dates.csv')
#colfunc=colorRampPalette(c("red","yellow","springgreen","royalblue"))
#cols=colfunc(100)
plot(dat$N~dat$Lon)

test=lmer(C~Lat+Lon+(1|SealID),data=dat) #add some metric of daily drift to this

uniqueseals=unique(dat$TOPPID)
#Color=brewer.pal(length(uniqueseals),"BrBG")
Color=c('#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd')
  
png(file='d15NoverTime.png',bg="white",units="in", width=5, height=5, res=600)
for(i in 1:length(uniqueseals)){
  dat_i=dat[dat$TOPPID==uniqueseals[i],]
  if(i==1){
    plot(dat_i$DaysSinceJune1,dat_i$N,
         xlab="Days Since June 01",ylab="d15N",
         ylim=c(12,16),xlim=c(0,300),type="l")
    par(new=TRUE)
    plot(dat_i$DaysSinceJune1,dat_i$N,
         xlab="",ylab="",xaxt="n",yaxt="n",
         ylim=c(12,16),xlim=c(0,300),type="p",pch=21,bg=Color[i])
  }else{
    par(new=TRUE)
    plot(dat_i$DaysSinceJune1,dat_i$N,
         xlab="",ylab="",xaxt="n",yaxt="n",
         ylim=c(12,16),xlim=c(0,300),type="l")
    par(new=TRUE)
    plot(dat_i$DaysSinceJune1,dat_i$N,
         xlab="",ylab="",xaxt="n",yaxt="n",
         ylim=c(12,16),xlim=c(0,300),type="p",pch=21,bg=Color[i])
  }
}
legend(x=220,y=15,legend=uniqueseals,pt.bg=Color,pch=21)
dev.off()

##########CARBON

png(file='d13CoverTime.png',bg="white",units="in", width=5, height=5, res=600)
for(i in 1:length(uniqueseals)){
  dat_i=dat[dat$TOPPID==uniqueseals[i],]
  if(i==1){
    plot(dat_i$DaysSinceJune1,dat_i$C,
         xlab="Days Since June 01",ylab="d13C",
         ylim=c(-18.5,-17),xlim=c(0,300),type="l")
    par(new=TRUE)
    plot(dat_i$DaysSinceJune1,dat_i$C,
         xlab="",ylab="",xaxt="n",yaxt="n",
         ylim=c(-18.5,-17),xlim=c(0,300),type="p",pch=21,bg=Color[i])
  }else{
    par(new=TRUE)
    plot(dat_i$DaysSinceJune1,dat_i$C,
         xlab="",ylab="",xaxt="n",yaxt="n",
         ylim=c(-18.5,-17),xlim=c(0,300),type="l")
    par(new=TRUE)
    plot(dat_i$DaysSinceJune1,dat_i$C,
         xlab="",ylab="",xaxt="n",yaxt="n",
         ylim=c(-18.5,-17),xlim=c(0,300),type="p",pch=21,bg=Color[i])
  }
}
legend(x=220,y=-17.5,legend=uniqueseals,pt.bg=Color,pch=21)
dev.off()

##########################
### ISOSCAPE CALCULATIONS SUMMARIZED
##########################

setwd("~/Publications/In Prep/NESE Isoscape Project/Isoscape Raw Data")
dat=read.csv('IsoscapeTripCalcs 2018_12_15.csv')
colfunc=colorRampPalette(c("red","yellow","springgreen","royalblue"))
cols=colfunc(100)
plot(dat$N~dat$Lon)





