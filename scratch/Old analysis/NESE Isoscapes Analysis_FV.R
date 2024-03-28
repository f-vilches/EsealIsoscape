setwd("C:/Doctorado/Research projects/Seals/FLO_NESE Isoscape Project/Isoscape Analysis 2021/Isoscape Analysis Dec")
dat=read.csv('NESEIsoscapeMetadata.csv')

LengthRecovery=dat$LengthRecovery 
LengthDeployment=dat$LengthDeployment
DeploymentDuration=dat$DeploymentDuration
Seal=dat$Seal

##### TO CALCULATE K FOR EACH WHISKER

ResultsMatrix=matrix(nrow=length(dat[,1]),ncol=3) #laying out an empty matrix, which will then be populated with the output of the for loop. For # of rows, the code uses a function to ask what is the length of the dataframe in column 1, and make that # the # of rows in the matrix.

for(i in 1:length(dat[,1])){ # i means whatever I want, it does not mean anything at this point, it is the variable that will take the value of the vector. The vector goes from 1 to the the last row of the dataframe.

#the 4 following lines are the variables of the matrix.
  
  asymptote=LengthDeployment[i] #assume that whiskers were at asymptotic length when plucked
  T=DeploymentDuration[i] #why T in red? R thinks is a shortcut for a logical vector (TRUE/FALSE), but it's not
  L=LengthRecovery[i]
  K=-(log(1-(L/asymptote))/T)
  
#the following 3 lines are the locations in the matrix for each variable. i means each iteration. E.g. for the 1st iteration: value of k will be placed in row i (which, for the 1st iteration is 1) and column 1; value of asymptote will be placed in row i (which, for the 1st iteration is 1), columnn 2. Row i for the 2nd itiration will be 2, and so on. 
  ResultsMatrix[i,1]=K
  ResultsMatrix[i,2]=asymptote
  ResultsMatrix[i,3]=T 
 
}

write.csv2(x = ResultsMatrix,file = "K_values.csv") #Export matrix

##### GRAPH GROWTH OF ALL WHISKERS

T0=0

png("WhiskerGrowth.png",width=5,height=5,res=600,units="in")
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