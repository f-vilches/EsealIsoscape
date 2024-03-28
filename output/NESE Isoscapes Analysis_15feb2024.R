#Script to match whisker stable isotope values with growth parameters and tracking coordinates----

#Project: Northern elephant seals Isoscapes
#Writer: Florencia Vilches
#Last modified: 17 Feb 2024

#read in R packages 
library("ggplot2")
library("ggfortify")	
library("lubridate")
library("mapdata") #package for pulling out map data
library("RColorBrewer")
library("reshape2")
library("xts")
#library("xtractomatic")
library(dplyr)
library(tidyverse)
require("httr")
require("ncdf4") 
require("sp")
library(lme4)
library(zoo)
library(here) 
library(ggformula)
library("lmerTest")
library("car")
library("brms")
library("bayesplot")
library("optimx")
library("ggthemes")
library("beepr")
library("egg")


Isotopes = read.csv("NESE_Isoscape_IsotopeData (R).csv", sep=";") # load csv file

#-----------------------------------------------------------------------------------#Chunk between dotted line is CHAOS!!!!!!!!! -Clean & streamline
  
#*Add new columns for and calculate each of the whisker growth model parameters----

Isotopes_WhiskerGrowth = Isotopes %>%
                             subset(SealID !="X822" & SealID != "5430"& SealID !="X 58"& SealID !="5055"& SealID !="B 66"& SealID !="9838"& SealID !="7019") %>% #remove non-nubbin seals
                             drop_na(d13C, LengthRecovery, RecoveryDate) %>% #remove rows with NA values in these columns
                             mutate(RecoveryDate= as.Date(RecoveryDate, format="%m/%d/%Y"),
                                    DeploymentDate= as.Date(DeploymentDate, format = "%m/%d/%Y"), 
                                    RecoveryDate= as.numeric(RecoveryDate), #it looks like a date, considered as a character. I think I cannot transform them directly from char to numeric if they look like dates.
                                    DeploymentDate= as.numeric(DeploymentDate), 
                                    Linf= 8.53 - 1.96*(WhiskerColumn-4)+2.06*(WhiskerRow-4)-0.56*(WhiskerRow-4)^2, #model to predict maximum whisker length (Beltran et al.2015)
                                    DeploymentDate= ifelse(Nubbin=="yes" |SealID== "7019"|SealID=="PM695",
                                                          DeploymentDate, 
                                                          (1/((exp(1))^-2.052821/(ifelse(Linf <= LengthRecovery,
                                                                                      LengthRecovery^0.9204057,
                                                                                      Linf^0.9204057))))*log(1-(ifelse(Linf <= LengthRecovery,
                                                                                                                        LengthRecovery/(LengthRecovery + 0.1), 
                                                                                                                        LengthRecovery/Linf)))
                                                                  +RecoveryDate),#For non-nubbins, the Deployment Date is solved for using Eq. 2 in Beltran et al. 2015, T0= 1/k*log(1-LengthRec/Linf)-RecoveryDate; plugging Lt and T as the recovery length and recovery date, respectively, Linf as above, and for K, Roxanne's equation to predict K from Linf: K <- exp(-2.052821) * asymp^(-0.9204057). Also, if the asymptotic length is shorter than or equal to the rec length, the log argument will be zero or negative. If so, replace the actual asymptotic length with the rec length plus 0.000001 (~0), to make the argument slightly positive.
                                                        k= -(log(1 - (ifelse(Linf <= LengthRecovery, 
                                                        LengthRecovery/(LengthRecovery + 0.1), 
                                                        LengthRecovery/Linf)))/(RecoveryDate-DeploymentDate)),#The argument for any log has to be a positive value. If the asymptotic length is shorter than or equal to the rec length, the argument will be zero or negative. If so, replace the actual asymptotic length with the rec length plus 0.000001 (~0), to make the argument slightly positive.
                                   DateWhenSegmentGrew = -1/k*log(1- (ifelse(Linf <= LengthRecovery, 
                                                                                     DistanceFromTip/LengthRecovery, 
                                                                                     DistanceFromTip/Linf)))                                                                                + DeploymentDate,#Age of the entire whisker for a given length (distance from tip), turned into date. The argument for any log has to be a positive value. If the asymptotic length is shorter than the rec length, the argument will be negative. If so, replace the actual asymptotic length with the rec length, to make the argument positive.
                                   DateWhenSegmentGrew= as.Date(DateWhenSegmentGrew),
                                   RecoveryDate= as.Date(RecoveryDate),
                                   DeploymentDate= as.Date(DeploymentDate), 
                                   DeploymentDuration= as.numeric(RecoveryDate - DeploymentDate), #from a time interval object (difftime) to a numeric object 
                                   SealDate=paste(SealID,DateWhenSegmentGrew,sep="-"),
                                   DateGrowthStart.forLauraTheseValuesWereListedUnderTime.= NULL,Comments= NULL, Weighter= NULL, Lat= NULL, Lon=NULL) #remove unnecessary columns 

Isotopes_WhiskerGrowth = Isotopes_WhiskerGrowth %>% #Seal 6564 is a special case; it's a nubbin but the length of the recovery whisker at the time of deployment was written on the datasheet (6.86 cm + base). The deployment length and date were plugged in the model to estimate T0 
  mutate(DeploymentDate = ifelse(SealID == "6564",
                                 1/((exp(1))^-2.052821/LengthRecovery^0.9204057)*log(1-(7.86/LengthRecovery))+17653,
                                 DeploymentDate),
         k = ifelse(SealID == "6564",
                    (exp(1))^-2.052821/LengthRecovery^0.9204057,
                    k),
         DateWhenSegmentGrew = ifelse(SealID == "6564",
                                      -1/k*log(1- (DistanceFromTip/LengthRecovery))+ DeploymentDate,
                                      DateWhenSegmentGrew),
          DeploymentDuration = ifelse(SealID == "6564",
                                     as.numeric(RecoveryDate - DeploymentDate),
                                     DeploymentDuration))

Isotopes_WhiskerGrowth= mutate(Isotopes_WhiskerGrowth, DeploymentDate=as.Date(DeploymentDate), DateWhenSegmentGrew=as.Date(DateWhenSegmentGrew))

#----------------------------------------------------------------------------------------------
#*For each whisker segment, calculate a mean d13C and d15N value from the previous 7 days (turnover rate)----

rolling_mean= function(x) 
{mean(x, na.rm = TRUE)}# Función para calcular la media de los últimos 7 días para cada fila

Isotopes_WhiskerGrowth= Isotopes_WhiskerGrowth %>% 
                            arrange("SealID", "DateWhenSegmentGrew") %>%
                            group_by(SealID) %>%
                            mutate(SevenDayAvg_d13C = sapply(seq_along(d13C), 
                                   function(i) {current_date= DateWhenSegmentGrew[i]
                                                last_seven_days= seq(current_date - 6,current_date, by = "day")
                                                rolling_mean(d13C[DateWhenSegmentGrew %in% last_seven_days])}),
                                   SevenDayAvg_d15N = sapply(seq_along(d15N), 
                                   function(i) {current_date= DateWhenSegmentGrew[i]
                                                last_seven_days= seq(current_date - 6,current_date, by = "day")
                                                rolling_mean(d15N[DateWhenSegmentGrew %in% last_seven_days])}))

#*Match DateWhenSegmentGrew with corresponding LAT and LON from MATLAB tracking data---- 

#Before running the code below, the Track_Best table in each diving data MATLAB file has to be exported as a csv file. Use MatLab code "export_TrackBest, DiveType, DiveStat_as_CSV_files" for that.

#Merge all TrackBest csv files of each seal into one dataframe: 

TrackBest_allseals= list.files(here("data/Tracking & diving"), pattern = "Track_Best", full.names = TRUE) %>% # produces a character vector of the names of files in the named directory or path, only for files that contain "Track_Best" in the file name. If TRUE, the directory path is prepended to the file names to give a relative file path. If FALSE, the file names (rather than paths) are returned.
                    map_df(~read_csv(.x) %>% # "map()" applies a function to each element of a list or atomic vector; "map_df()" returns a data frame; argument ".x" is a list or atomic vector
                    mutate(SealID= as.character(SealID), #to prevent the error "Can't Combine X <character> and X <double>" for SealID
                           Lon= Long)) #Longitude is abbreviated as "Lon" everywhere but for some reason it shows as "Long" in the TrackBest tables exported as CSV, this line fixes it

TrackBest_allseals$DateUTC=as.Date(TrackBest_allseals$DateUTC, format="%d-%b-%Y") #dates from character objects to date objects and remove the H-M-S component of the date-time; Date UTC column is present in some TrackBest tables, others have only JulDate. For the latter, this line of code is not necessary

# Create the function that will convert MATLAB Julian date to R date:

matlab2POS = function(x, timez = "UTC") {
  days = x - 719529 	# 719529 = days from 1-1-0000 to 1-1-1970
  secs = days * 86400 # 86400 seconds in a day
  
  return(as.POSIXct(strftime(as.POSIXct(secs, origin = '1970-1-1', 
                                        tz = 'UTC'), format = '%Y-%m-%d %H:%M', 
                             tz = 'UTC', usetz = FALSE), tz = timez))
}

#Reduce the dataframe to 1-hit per day (median of date/time):

TrackBest_allseals_subsampled= TrackBest_allseals %>% 
                               mutate(DateWhenSegmentGrew=as.Date(matlab2POS(JulDate)),#Convert the Julian dates from MATLAB format into R format, using the function above
                                      SealDate=paste(SealID,DateWhenSegmentGrew,sep="-")) %>% 
                               group_by(SealDate) %>% 
                               summarize(Lat=median(Lat), Lon=median(Lon)) 

#For each DateWhenSegmentGrew (in df Isotopes), assign the corresponding latitude and lognitude (in df TrackBest_allseals):

IsotopesLatLon=Isotopes_WhiskerGrowth %>%
                        left_join(TrackBest_allseals_subsampled, by= "SealDate")%>% #left_join function retains only rows for which the first (left, x, Isotopes) data frame has complete data; the second (right, y) data frame is TrackBest
                        arrange("SealID", "Segment") %>% #After merging, the rows shuffle. Sort them back by Seal ID and Segment.
                        mutate(DaysSinceDeployment= as.numeric(DateWhenSegmentGrew - DeploymentDate),#Will be used as the time reference on the x-axis. Day of year is useless because Jan 1st is at the end of the deployment and whisker growth. Transform it from a time interval object (difftime) to a numeric object
                               Lon360= ifelse(Lon>0, 
                                              Lon-360,
                                              Lon)) 

#*Match DateWhenSegmentGrew with corresponding DRIFT RATE from MATLAB tracking data----

#Merge all DiveType csv files of each seal into one dataframe - To get the drift rate: 

DiveType_allseals= list.files(here("data/Tracking & diving"), pattern = 'DiveType',full.names = TRUE) %>%  # "list.files()" produces a character vector of the names of files in the named directory or path, only for files that contain "DiveType" in the file name. If TRUE, the directory path is prepended to the file names to give a relative file path. If FALSE, the file names (rather than paths) are returned.
                  map_df(~read_csv(.x) %>% # "map()" applies a function to each element of a list or atomic vector; "map_df()" returns a data frame; argument ".x" is a list or atomic vector
                  mutate(SealID = as.character(SealID))) #to prevent the error "Can't Combine X <character> and X <double>" for SealID

#Merge all DiveStat csv files of each seal into one dataframe - To get the date for each drift rate: 

DiveStat_allseals= list.files(here("data/Tracking & diving"), pattern = 'DiveStat',full.names = TRUE) %>%  # "list.files()" produces a character vector of the names of files in the named directory or path, only for files that contain "DiveStat" in the file name. 
                   map_df(~read_csv(.x) %>% # "map()" applies a function to each element of a list or atomic vector; "map_df()" returns a data frame; argument ".x" is a list or atomic vector
                   mutate(SealID = as.character(SealID))) #to prevent the error "Can't Combine X <character> and X <double>" for SealID

#Merge DiveType & DiveStat:

DriftRate=DiveType_allseals %>%
          select(SealID, DiveNumber, DriftRate) %>% #columns to keep from DiveType
          left_join(select(DiveStat_allseals, SealID, DiveNumber, JulDate), by=c("SealID", "DiveNumber")) %>% #columns to keep from DiveStat; left join retains only rows for which the first (left, x, DiveType) data frame has complete data; the second (right, y) data frame is DiveStat
          mutate(DateWhenSegmentGrew=as.Date(matlab2POS(JulDate), format = '%Y-%m-%d')) %>% # Convert the Julian dates from MATLAB format into R format, using the matlab2POS function created above
          filter(abs(DriftRate - mean(DriftRate)) <= 3 * sd(DriftRate)) %>% #eliminate all data points that are more than >3 SD from the mean
          group_by(SealID,DateWhenSegmentGrew) %>% 
          summarize(DriftRate=mean(DriftRate)) %>%
          mutate(FiveDayDriftRate = rollapply(DriftRate, width = 5, FUN = mean, align = "center", fill = NA)) #Reduce the dataframe to a mean drift rate per day and then to one drift rate every 5 days by calculating a 3-day centered moving average for each seal

#For each DateWhenSegmentGrew (in df IsotopesWithLatLon), assign the corresponding drift rate (in df DriftRate):

IsotopesLatLonDrift = IsotopesLatLon %>%
  mutate(datekey = as.character(DateWhenSegmentGrew)) %>% #turn Date into a character and put it in a new column bc the actual date had hidden seconds and was not properly matched
  left_join(transmute( 
    DriftRate, 
    SealID, 
    datekey = as.character(DateWhenSegmentGrew), 
    FiveDayDriftRate),  #columns to keep from DriftRate; left join retains only rows for which the first (left, x, IsotopesLatLon) data frame has complete data; the second (right, y) data frame is DriftRate 
    by = c("SealID", "datekey"
    )) %>% 
  select(-datekey)
                                       
#*Match DateWhenSegmentGrew with corresponding REPRODUCTIVE STATE----

ReproductiveState = read.csv("NESE_Reproductive state.csv", sep=";") # load csv file
IsotopesLatLonDriftRepr=merge(IsotopesLatLonDrift,
                                       ReproductiveState) #Add reproductive state data

#*Exploratory plots----

#**Histograms of the isotope values (response variables distribution)----

par(mfrow = c(1, 2), mai = c(0.9, 0.8, 0.2, 0.8))  # Set the number of rows and columns in the plotting device & adjust the margin to create more space
hist(IsotopeLatLonDriftRepr$d13C, main="Carbon", xlab="d13C (‰)")
hist(IsotopeLatLonDriftRepr$d15N, main="Nitrogen", xlab ="d15N (‰)")

#**Drift Rate Vs Days since deployment to check for outliers

DriftRate= DriftRate %>% 
           group_by(SealID) %>%
           mutate(DeploymentDuration= row_number()) #For each seal, turn each row into day of deployment

ggplot(data = DriftRate, aes(x = DeploymentDuration, y = DriftRate)) +
  geom_point(size=0.2) +
  geom_smooth(size=0.7, color="red") +
  labs(x="Deployment Duration (days)", y="Daily Drift Rate (m/s)") +
  facet_wrap(~ SealID, scales = "free") +
  theme_classic() 

ggplot(data = DriftRate, aes(x = DeploymentDuration, y = ThreeDayCenteredAvgDriftRate)) +
  geom_point(size=0.2) +
  geom_smooth(size=0.7,color="red") +
  labs(x="Deployment Duration (days)", y="3-Day Mov Avg Drift Rate (m/s)") +
  facet_wrap(~ SealID, scales = "free") +
  theme_classic() 

#**Relationships between responses and predictors----

#Plot isotope values versus all numerical predictors
pairs(d13C ~ Lat+Lon360+FiveDayDriftRate,data=IsotopesLatLonDriftRepr,panel=panel.smooth)
pairs(d15N ~ Lat+Lon360+FiveDayDriftRate,data=IsotopesLatLonDriftRepr,panel=panel.smooth)

#Plot isotope values versus time since deployment

ggplot(data=IsotopesLatLonDriftRepr, aes(x=DaysSinceDeployment, y=d13C)) +
  geom_line(size=.8) +
  facet_wrap(~ SealID) +
  theme_classic()

ggplot(data=IsotopesLatLonDriftRepr, aes(x=DaysSinceDeployment, y=d15N)) +
  geom_line(size=.8) +
  facet_wrap(~ SealID) +
  theme_classic()

#Plot isotope values versus lat, lon, drift by seal

ggplot(data=IsotopesLatLonDriftRepr, aes(x=Lat, y=d13C, color=SealID)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_few()

ggplot(data=IsotopesLatLonDriftRepr, aes(x=FiveDayDriftRate, y=d13C, color=SealID)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_few()

ggplot(data=IsotopesLatLonDriftRepr, aes(x=Lon360, y=d15N, color=SealID)) +
  geom_point() +
  stat_smooth(method="lm") +
  facet_wrap(~SealID) +
  theme_few() #Para lon, casi todas las seals tienen pendiente +

ggplot(data=IsotopesLatLonDriftRepr, aes(x=FiveDayDriftRate, y=d15N, color=SealID)) +
  geom_point() +
  stat_smooth(method="lm") +
  facet_wrap(~SealID) +
  theme_few() #Para DR, ~ proporciones de +, - y nula

#Plot isotope values versus reproductive state

ggplot(data=IsotopesLatLonDriftRepr, aes(x=ReprState, y=d13C))+
  geom_boxplot()+ 
  geom_jitter(height=0, color='black', fill='grey50', shape=21, alpha= 0.2)+ #alpha changes the transparency
  theme_classic()

ggplot(data=IsotopesLatLonDriftRepr, aes(x=ReprState, y=d15N))+
  geom_boxplot()+ 
  geom_jitter(height=0, color='black', fill='grey50', shape=21, alpha= 0.2) +
  theme_classic()

# **Isotope maps----

yrange=c(33,65) #lat range of all the seals
xrange=c(-183,-110) #lon range of all the seals

w=map_data("world",ylim=yrange,xlim=xrange) #world map data

states=map_data("state") #states data

# ***For CARBON----

assign_range = function(d13C) { #function to group isotope values into ranges to better visualize them on the map
  range_limits = c(-19, -18.5, -18, -17.5, -17)
  
  if (d13C >= range_limits[1] && d13C <= range_limits[2] ) { #assign ranges based on the limits
    return("Low")
  } else if (d13C > range_limits[2] && d13C <= range_limits[3]) {
    return("Low-medium")
  } else if (d13C > range_limits[3] && d13C <= range_limits[4]) {
    return("Medium")
  } else if (d13C > range_limits[4] && d13C <= range_limits[5]) {
    return("Medium-high")
  } else {
    return("High")
  }
}

IsotopesLatLon$d13C_ranges = sapply(IsotopesLatLon$d13C, assign_range) #apply the range assignment function to the column with d13C values

ggplot()+labs(y= "Latitude (deg)", x = "lonitude (deg)")+
  geom_polygon(data=w,aes(x=lon,y=lat,group=group),fill="grey30")+
  coord_fixed(xlim=xrange,ylim=yrange)+ #Cartesian coordinate system with a fixed aspect ratio.
  geom_point(data=IsotopesLatLon,shape=21, size=2, color="grey30", aes(x=lon,y=Lat,fill=d13C_ranges)) +
  scale_fill_manual(name="δ13C",
                       limits = c("Low", "Low-medium", "Medium", "Medium-high", "High"), #give the items in the legend the order I want
                       labels = c("-19 to-18.5", "-18.4 to -18", "-17.9 to -17.5", "-17.4 to -17", "-16.9 to -16"), #match labels to the corresponding item on the legend 
                       values = c("#feeada","#fca562","#f96e03","#984302","#381900"))+ #change the color of the datapoints
  theme(legend.position = c(0.9, 0.76)) #position the legend inside the plotting area; x and y are the coordinates of the legend box, should be between 0 and 1
  
# ***For NITROGEN----

assign_range = function(d15N) { #function to group isotope values into ranges to better visualize them on the map
  range_limits = c(11, 12.5, 14, 15.5, 17)
  
  if (d15N >= range_limits[1] && d15N < range_limits[2] ) { #assign ranges based on the limits
    return("Low")
  } else if (d15N >= range_limits[2] && d15N < range_limits[3]) {
    return("Medium")
  } else if (d15N >= range_limits[3] && d15N < range_limits[4]) {
    return("Medium-high")
   } else {
    return("High")
  }
}

IsotopesLatLon$d15N_ranges = sapply(IsotopesLatLon$d15N, assign_range) #apply the range assignment function to the column with d13C values
  
ggplot()+labs(y= "Latitude (deg)", x = "lonitude (deg)")+
  geom_polygon(data=w,aes(x=lon,y=lat,group=group),fill="grey30")+
  coord_fixed(xlim=xrange,ylim=yrange)+ #Cartesian coordinate system with a fixed aspect ratio.
  geom_point(data=IsotopesLatLon,shape=21, size=2, color="grey30", aes(x=lon,y=Lat,fill=d15N_ranges))+
  scale_fill_manual(name="δ15N",
                  limits = c("Low", "Medium","Medium-high", "High"), #give the items in the legend the order I want
                  labels = c("11 to 12.4", "12.5 to 13.9", "14 to 15.4", "15.5 to 17.3"), #match labels to the corresponding item on the legend 
                  values = c("#dff1e1","#6cc173","#2e6c33","#102812"))+ #change the color of the datapoints
  theme(legend.position = c(0.9, 0.76)) #position the legend inside the plotting area; x and y are the coordinates of the legend box, should be between 0 and 1


#*Fit a Linear mixed effect model----

IsotopesLatLonDriftRepr = IsotopesLatLonDriftRepr[complete.cases(IsotopesLatLonDriftRepr$FiveDayDriftRate, IsotopesLatLonDriftRepr$Lon360, IsotopesLatLonDriftRepr$ReprState),] #remove any rows with missing values in the FiveDayDriftRate, Lon360, and ReprState columns, which ensures that all the models below are fitted to the same subset of the data

#**For CARBON----

#*Fit the model:

LMM_full_Carbon = lmer(d13C~Lat + Lon360 + FiveDayDriftRate + ReprState + (1| SealID), 
                  data = IsotopesLatLonDriftRepr); summary(LMM_full_Carbon)

#Test assumptions:

shapiro.test(resid(LMM_full_Carbon)) #not normal
plot(LMM_full_Carbon) #plot the residuals
qqnorm(residuals(LMM_full_Carbon)) 
qqline(residuals(LMM_full_Carbon)) 

#Normal QQ plot, has a heavy tail, so a linear model with a brm function and a student distribution is more appropiate:

LMM_full_Carbon_brm2 = brm(d13C~Lat + Lon360 + FiveDayDriftRate + ReprState + (Lat| SealID),
                           data = IsotopesLatLonDriftRepr,family = student); summary(LMM_full_Carbon_brm2);beep(3); fixef(LMM_full_Carbon_brm2) #prints more digits for easier viewing

#Assess model's predicting abilities:

bayesplot_grid(pp_check(
  LMM_full_Carbon_brm2, type= "stat", stat="mean")) #observed mean (black) relative to the means from simulated data (blue) plugged into the model
  
bayesplot_grid(pp_check(
  LMM_full_Carbon_brm2, ndraws=50)) #density of observed responses (black) relative to the density of responses from simulated data (blue) plugged into the model

#**For NITROGEN----

#Rescale predictors to get rid of the very large eigenvalue: 

IsotopesLatLonDriftRepr= IsotopesLatLonDriftRepr %>%  
  mutate(LatRescaled= scale(Lat),
         LonRescaled= scale(Lon360)) 

#Fit the model:

GLMM_full_Nitro = glmer(d15N~LatRescaled + LonRescaled + FiveDayDriftRate + ReprState  + (1| SealID),
                        family = Gamma (link='log'),
                        data = IsotopesLatLonDriftRepr); 
                        summary(GLMM_full_Nitro)

#Extract the fixed effects coefficients (slopes) and convert them from the log link scale to the response scale:

fixed_effects= fixef(GLMM_full_Nitro); slopes_response_scale= exp(fixed_effects); print(slopes_response_scale) 

#Create “new data” data frame
nd2=data.frame(LatRescaled=rep(-3:4, times=23), LonRescaled=rep(-3:4, times=23),FiveDayDriftRate=rep(-3:4, times=23), ReprState=rep(c("Pregnant","Non-pregnant"),each=184), SealID="3190")

#Add predicted values
nd2$pred_nitro=predict(GLMM_full_Nitro,newdata=nd2,type="response", re.form=NA)

#Generate confidence intervals: 
nd2$Nitro_link=predict(GLMM_full_Nitro,newdata=nd2,re.form=NA,type="link") #predict mean values on link/log scale
pf1 = function(fit) {   predict(fit, nd2) } #function for bootstrapping
bb=bootMer(GLMM_full_Nitro,nsim=50,FUN=pf1,seed=69) #bootstrap to estimate uncertainty in predictions
nd2$SE= apply(bb$t, 2, sd) #Calculate Ses from bootstrap samples on link scale
nd2$d15N= exp(nd2$Nitro_link) #predicted mean values on response scale
nd2$pSE= exp(nd2$Nitro_link+nd2$SE) #predicted mean + 1 SE on response scale
nd2$mSE= exp(nd2$Nitro_link-nd2$SE) # predicted mean - 1 SE on response scale

ggplot(data=IsotopesLatLonDriftRepr,aes(x=LonRescaled,y=d15N,color=SealID)) +
         geom_point(size=3)+ 
         geom_ribbon(data=nd2,aes(x=LonRescaled,ymin=mSE,ymax=pSE,color=SealID), alpha=0.1, linetype=0)+
         theme_few()+
         geom_line(data=nd2,aes(x=LonRescaled,y=d15N,color=SealID))+
         labs(x = "LonRescaled",y="δ15N")

#Fit a second model, with DR as a random effect:

GLMM_full_Nitro2 = glmer(d15N~LatRescaled + LonRescaled + FiveDayDriftRate + ReprState  + (FiveDayDriftRate|SealID),
                        family = Gamma (link='log'),
                        data = IsotopesLatLonDriftRepr); 
                        summary(GLMM_full_Nitro2)

#Extract the fixed effects coefficients (slopes) and convert them from the log link scale to the response scale: 
fixed_effects= fixef(GLMM_full_Nitro2); slopes_response_scale= exp(fixed_effects); print(slopes_response_scale) 

#*Actual figures----

#**For CARBON----

CarbonVsLat= ggplot(data=IsotopesLatLonDriftRepr, aes(x=Lat, y=d13C, color=SealID)) +
  geom_point(size= 1, show.legend = FALSE) +
  stat_smooth(method="lm", alpha = 0.07, linewidth= 0.6, show.legend = FALSE, aes(fill=SealID)) + # Change the color of the confidence intervals to blue) + #plot the model and confidence intervals
  labs(x="Latitude (degrees)", y=bquote("Carbon isotope values, δ" ^13~"C"~"(‰)")) +
  theme_few(); CarbonVsLat

CarbonVsDrift= ggplot(data=IsotopesLatLonDriftRepr, aes(x=FiveDayDriftRate, y=d13C, color=SealID)) +
  geom_point(size= 1, show.legend = FALSE) +
  stat_smooth(method="lm", alpha = 0.06, linewidth= 0.6, show.legend = FALSE, aes(fill=SealID)) + # Change the color of the confidence intervals to blue) + #plot the model and confidence intervals
   labs(x="Drift rate (m/s)", y=bquote("Carbon isotope values, δ" ^13~"C"~"(‰)")) +
  theme_few(); CarbonVsDrift

IsotopeRelation= ggarrange(CarbonVsLat, CarbonVsDrift, 
                           ncol = 2, nrow = 1) #Combine carbon plots

#**For NITROGEN----



#*Test assumptions----IGNORE

# Carbon Vs Lat
lmCarbonLat= lm(d13C~Lat, data=IsotopesLatLon);summary(lmCarbonLat) #Fit the model
shapiro.test(resid(lmCarbonLat)) #test normal distribution of residuals; normal distribution is not met
lmtest::bptest(lmCarbonLat) # test constant variance of residuals; constant variance is met
dwtest(lmCarbonLat) # formal test for independence of residuals (Durbin Watson Test; H0: autocorrelation does not exist, H1: autocorrelation exists) 
pacf(resid(lmCarbonLat)) # graph test for independence of residuals  - autocorrelation exists