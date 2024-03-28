setwd("C:/Users/roxan/Desktop")
dat=read.csv('LauraIsoscapeMetadata.csv')

LengthRecovery=dat$LengthRecovery
LengthDeployment=LengthRecovery
DeploymentDuration=dat$DeploymentDuration
Seal=dat$Seal

##### TO CALCULATE K FOR EACH WHISKER

ResultsMatrix=matrix(nrow=length(dat[,1]),ncol=3)

for(i in 1:length(dat[,1])){ #for each whisker
  asymptote=LengthDeployment[i] #assume that whiskers were at asymptotic length when plucked
  T=DeploymentDuration[i] 
  L=LengthRecovery[i]
  K=-(log(1-(L/asymptote))/T)
  
  ResultsMatrix[i,1]=K
  ResultsMatrix[i,2]=asymptote
  ResultsMatrix[i,3]=T
  
}

##### GRAPH GROWTH OF ALL WHISKERS

T0=0

png("LauraWhiskerGrowth.png",width=5,height=5,res=600,units="in")
par(mfrow=c(1,1),xpd=FALSE,mar=c(5,4,1,1))
for(i in 1:length(dat[,1])){
  k=0.0132
  asymp=dat$LengthRecovery[i]
  x_temp=seq(1,dat$DeploymentDuration[i],.1)
  y_temp=asymp*(1-exp(-k*(x_temp-T0)))
  if(i==1){
    plot(x_temp,y_temp,type="l",lwd=.5,xlab="Time (days)",ylab="Length (cm)",xlim=c(0,365),ylim=c(0,14))
    abline(v=dat$DeploymentDuration,col="darkgrey")
  }else{
    par(new=TRUE)
    plot(x_temp,y_temp,type="l",lwd=.5,xlab="",ylab="",xlim=c(0,365),ylim=c(0,14),xaxt="n",yaxt="n")
  }
}

dev.off()


### TO PLOT CARBON OR NITROGEN ACROSS SEGMENTS
for(i in 1:length(AllSeals)){
  jpeg(file = paste(AllSeals[i],"_Nitrogen",".jpeg", sep=""))
  
  x_temp=dat$Segment[dat$Seal==AllSeals[i]]
  y_temp=dat$d15N....[dat$Seal==AllSeals[i]]
  plot(x_temp,y_temp,ylim=c(12,17),xlim=c(24,0),xlab="Segment #",ylab="d15N",col="black",type="b",main=AllSeals[i])
  dev.off()
}