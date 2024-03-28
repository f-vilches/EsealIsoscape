setwd("~/Collaborations/Laura Christianson/Isoscape 2019_10_31")
dat=read.csv('LauraIsoscapeMetadata.csv')
dat2=read.csv('Isotopes_FillIn.csv')

TOPPIDs=unique(dat$TOPPID)

K=0.0132

for(i in 1:nrow(dat2)){
T0=dat$T0[dat$TOPPID==dat2$TOPPID[i]]
asymp=dat$LengthRecovery[dat$TOPPID==dat2$TOPPID[i]]
Lt=asymp-dat2$CmFromBase[i]#distance from tip?
#y=asymp*(1-exp(-K*(x-T0))) #to calculate length as a function of date
#to calculate date as a function of length
dat2$T[i]=-1/K*log(1-Lt/asymp)+T0

#check with spreadsheet
}
