#Script to match whisker stable isotope values with growth parameters and tracking coordinates----

#Project: Northern elephant seals Isoscapes
#Writer: Florencia Vilches
#Last modified: 6 Feb 2024

#read in R packages 
library("ggplot2")
library("ggfortify")	
library("lubridate")
library("mapdata")
library("RColorBrewer")
library("reshape2")
library("xts")
#library("xtractomatic")
library(dplyr)
library(tidyverse)
require("httr")
require("ncdf4") 
require("sp")
library(mapdata) #package for pulling out map data
library(lme4)
library(zoo)
library(here) #enables easy file referencing by using the top-level directory of a file project to easily build file paths. In contrast to setwd(), which is fragile and dependent on the way you order your files on your computer.
library(ggformula)

Isotopes = read.csv("NESE_Isoscape_IsotopeData (R).csv", sep=";") # load csv file

-----------------------------------------------------------------------------------#Chunk between dotted line is CHAOS!!!!!!!!! -Clean & streamline
  
#*Add new columns for and calculate each of the whisker growth model parameters----

Isotopes_WhiskerGrowth = Isotopes %>%
                             subset(SealID !="X822" & SealID != "5430"& SealID !="X 58"& SealID !="5055"& SealID !="B 66"& SealID !="9838") %>% #remove non-nubbin seals
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

----------------------------------------------------------------------------------------------
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
                               mutate(DateWhenSegmentGrew=as.Date(matlab2POS(JulDate)), #Convert the Julian dates from MATLAB format into R format, using the function above
                                      SealDate=paste(SealID,DateWhenSegmentGrew,sep="-")) %>% 
                               group_by(SealDate) %>% 
                               summarize(Lat=median(Lat), Lon=median(Lon)) 

#For each DateWhenSegmentGrew (in df Isotopes), assign the corresponding latitude and lognitude (in df TrackBest_allseals):

IsotopesWithLatLon=Isotopes_WhiskerGrowth %>%
                        left_join(TrackBest_allseals_subsampled, by= "SealDate")%>% #left_join function retains only rows for which the first (left, x, Isotopes) data frame has complete data; the second (right, y) data frame is TrackBest
                        arrange("SealID", "Segment") %>% #After merging, the rows shuffle. Sort them back by Seal ID and Segment.
                        mutate(DaysSinceDeployment= as.numeric(DateWhenSegmentGrew - DeploymentDate)) #Will be used as the time reference on the x-axis. Day of year is useless because Jan 1st is at the end of the deployment and whisker growth. Transform it from a time interval object (difftime) to a numeric object

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
          group_by(SealID,DateWhenSegmentGrew) %>% 
          summarize(DriftRate=mean(DriftRate)) %>%
          mutate(ThreeDayCenteredAvgDriftRate = rollapply(DriftRate, width = 3, FUN = mean, align = "center", fill = NA)) #Reduce the dataframe to a mean drift rate per day and then to one drift rate every 3 days by calculating a 3-day centered moving average for each seal

#For each DateWhenSegmentGrew (in df IsotopesWithLatLon), assign the corresponding drift rate (in df DriftRate):

IsotopesWithLatLonDrift = IsotopesWithLatLon %>%
  mutate(datekey = as.character(DateWhenSegmentGrew)) %>% #turn Date into a character and put it in a new column bc the actual date had hidden seconds and was not properly matched
  left_join(transmute( 
    DriftRate, 
    SealID, 
    datekey = as.character(DateWhenSegmentGrew), 
    ThreeDayCenteredAvgDriftRate),  #columns to keep from DriftRate; left join retains only rows for which the first (left, x, IsotopesLatLon) data frame has complete data; the second (right, y) data frame is DriftRate 
    by = c("SealID", "datekey"
    )) %>% 
  select(-datekey)
                                       
#*Match DateWhenSegmentGrew with corresponding REPRODUCTIVE STATE----

ReproductiveState = read.csv("NESE_Reproductive state.csv", sep=";") # load csv file
IsotopeWithLatLonDriftRepr=merge(IsotopeWithLatLonDrift,
                                       ReproductiveState) #Add reproductive state data

#*Exploratory plots----

#**Histograms of the isotope values (response variables distribution)----

par(mfrow = c(1, 2), mai = c(0.9, 0.8, 0.2, 0.8))  # Set the number of rows and columns in the plotting device & adjust the margin to create more space
hist(IsotopeWithLatLonDriftRepr$d13C, main="Carbon", xlab="d13C (‰)")
hist(IsotopeWithLatLonDriftRepr$d15N, main="Nitrogen", xlab ="d15N (‰)")

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
pairs(d13C ~ Lat+Lon+DriftRate,data=IsotopesWithLatLonDriftRepr,panel=panel.smooth)
pairs(d15N ~ Lat+Lon+DriftRate,data=IsotopesWithLatLonDriftRepr,panel=panel.smooth)

#Plot isotope values versus time since deployment

ggplot(data=IsotopesWithLatLon, aes(x=DaysSinceDeployment, y=d13C)) +
  geom_line(size=.8) +
  theme_classic()

ggplot(data=IsotopesWithLatLon, aes(x=DaysSinceDeployment, y=d15N)) +
  geom_line(size=.8) +
  facet_wrap(~ SealID) +
  theme_classic()

#Plot isotope values versus lat and lon by seal

ggplot(data=IsotopesWithLatLon, aes(x=Lat, y=d13C, color=SealID)) +
  geom_line(size=.8) +
  theme_classic()

ggplot(data=IsotopesWithLatLon, aes(x=Lat, y=d15N, color=SealID)) +
  geom_line(size=.8) +
  theme_classic()

#Plot isotope values versus reproductive state

ggplot(data=IsotopesWithLatLon, aes(x=ReprState, y=d13C))+
  geom_boxplot()+ 
  geom_jitter(height=0, color='black', fill='grey50', shape=21, alpha= 0.2)+ #alpha changes the transparency
  theme_classic()

ggplot(data=IsotopesWithLatLon, aes(x=ReprState, y=d15N))+
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

IsotopesWithLatLon$d13C_ranges = sapply(IsotopesWithLatLon$d13C, assign_range) #apply the range assignment function to the column with d13C values

ggplot()+labs(y= "Latitude (deg)", x = "lonitude (deg)")+
  geom_polygon(data=w,aes(x=lon,y=lat,group=group),fill="grey30")+
  coord_fixed(xlim=xrange,ylim=yrange)+ #Cartesian coordinate system with a fixed aspect ratio.
  geom_point(data=IsotopesWithLatLon,shape=21, size=2, color="grey30", aes(x=lon,y=Lat,fill=d13C_ranges)) +
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

IsotopesWithLatLon$d15N_ranges = sapply(IsotopesWithLatLon$d15N, assign_range) #apply the range assignment function to the column with d13C values
  
ggplot()+labs(y= "Latitude (deg)", x = "lonitude (deg)")+
  geom_polygon(data=w,aes(x=lon,y=lat,group=group),fill="grey30")+
  coord_fixed(xlim=xrange,ylim=yrange)+ #Cartesian coordinate system with a fixed aspect ratio.
  geom_point(data=IsotopesWithLatLon,shape=21, size=2, color="grey30", aes(x=lon,y=Lat,fill=d15N_ranges))+
  scale_fill_manual(name="δ15N",
                  limits = c("Low", "Medium","Medium-high", "High"), #give the items in the legend the order I want
                  labels = c("11 to 12.4", "12.5 to 13.9", "14 to 15.4", "15.5 to 17.3"), #match labels to the corresponding item on the legend 
                  values = c("#dff1e1","#6cc173","#2e6c33","#102812"))+ #change the color of the datapoints
  theme(legend.position = c(0.9, 0.76)) #position the legend inside the plotting area; x and y are the coordinates of the legend box, should be between 0 and 1


#*Fit a generalized linear mixed effect model----

GLMM_full = glmer(d13C~Lat + Lon + DriftRate + ReprState + 
                  DriftRate*Lat+ DriftRate*Lon + DriftRate*ReprState + 
                  (1| SealID), 
                  data = IsotopesWithLatLonDriftRepr)

GLMM_full = glmer(d15N~Lat + Lon + DriftRate + ReprState + 
                  DriftRate*Lat+ DriftRate*Lon + DriftRate*ReprState + 
                  (1| SealID),
                  family = Gamma,
                  data = IsotopesWithLatLonDriftRepr)
#"calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated; please call lmer() directly - ChatGPT: using glmer with family = gaussian (identity link) is deprecated for fitting linear mixed-effects models and advises using lmer directly for this purpose. This change was introduced to make the syntax more explicit and avoid confusion between linear and generalized linear mixed-effects models. To address this, you should use the lmer function instead of glmer when fitting a linear mixed-effects model with a Gaussian distribution (identity link). Here's how you can modify your code

#*Test assumptions----

# Carbon Vs Lat
lmCarbonLat= lm(d13C~Lat, data=IsotopesWithLatLon);summary(lmCarbonLat) #Fit the model
shapiro.test(resid(lmCarbonLat)) #test normal distribution of residuals; normal distribution is not met
lmtest::bptest(lmCarbonLat) # test constant variance of residuals; constant variance is met
dwtest(lmCarbonLat) # formal test for independence of residuals (Durbin Watson Test; H0: autocorrelation does not exist, H1: autocorrelation exists) 
pacf(resid(lmCarbonLat)) # graph test for independence of residuals  - autocorrelation exists

#Carbon Vs lon
lmCarbonlon= lm(d13C~lon, data=IsotopesWithLatLon);summary(lmCarbonlon) #Fit the model
shapiro.test(resid(lmCarbonlon)) #test normal distribution of residuals; normal distribution is not met
lmtest::bptest(lmCarbonlon) # test constant variance of residuals; constant variance is met
dwtest(lmCarbonlon) # formal test for independence of residuals (Durbin Watson Test; H0: autocorrelation does not exist, H1: autocorrelation exists) 
pacf(resid(lmCarbonlon)) # graph test for independence of residuals  - autocorrelation exists

#Nitro Vs Lat
lmNitroLat= lm(d15N~Lat, data=IsotopesWithLatLon);summary(lmNitroLat) #Fit the model
shapiro.test(resid(lmNitroLat)) #test normal distribution of residuals; normal distribution is not met
lmtest::bptest(lmNitroLat) # test constant variance of residuals; constant variance is met
dwtest(lmNitroLat) # formal test for independence of residuals (Durbin Watson Test; H0: autocorrelation does not exist, H1: autocorrelation exists) 
pacf(resid(lmNitroLat)) # graph test for independence of residuals  - autocorrelation exists

#Nitro Vs lon
lmNitrolon= lm(d15N~lon, data=IsotopesWithLatLon);summary(lmNitrolon)  #Fit the model
shapiro.test(resid(lmNitrolon)) #test normal distribution of residuals; normal distribution is not met
lmtest::bptest(lmNitrolon) # test constant variance of residuals; constant variance is not met
dwtest(lmNitrolon) # formal test for independence of residuals (Durbin Watson Test; H0: autocorrelation does not exist, H1: autocorrelation exists) 
pacf(resid(lmNitrolon)) # graph test for independence of residuals  - autocorrelation exists