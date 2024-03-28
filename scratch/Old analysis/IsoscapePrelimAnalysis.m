format compact
format long
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Read in data, plot, setup for Iknos %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For all animals, define limits of plot
Lon1=-160
Lon2=-110
Lat1=30
Lat2=60

%For main map figure
figure;  %Open new figure
hold on %Keep the plot open so data can be added
axis([Lon1 Lon2 Lat1 Lat2]) %Change axis limits

filename='NESE Isoscape Unfiltered Tracks.csv';
[PTT,Month,Day,Year,Hour,Minute,Latitude,Longitude]=textread(filename,'%f , %f / %f / %f %f : %f , %f , %f','headerlines',1);

%Convert track date to julian date
TrackJulDate=datenum(Year,Month,Day);
    
Seals=unique(PTT);

%Load in isotope data with dates
filename='NESE Isoscape Isotope Data Dates.csv';
[PTT_Iso,Carbon,Nitrogen,Month_Iso,Day_Iso,Year_Iso]=textread(filename,'%f ,%f, %f, %f / %f / %f ','headerlines',1);
IsoJulDate=datenum(Year_Iso,Month_Iso,Day_Iso);

for i = 1:max(size(Seals))  %For each animal
    
    %Subset data for this seal
    Lon=Longitude(PTT==Seals(i));
    Lat=Latitude(PTT==Seals(i));
    Date=TrackJulDate(PTT==Seals(i));
     
    DateIso=IsoJulDate(PTT_Iso==Seals(i));  %Find number of isotope points for this seal 
    Carb=Carbon(PTT_Iso==Seals(i));
    Nitr=Nitrogen(PTT_Iso==Seals(i));
    
    for j=1:max(size(DateIso)) %For each isotope sample
   
    %Find latitude/longitude for each of the isotope points

    tmp=abs(DateIso(j)-Date); %Find row number of closest time
    [idx idx] = min(tmp); %index of closest value

    %TimeClosest=Date(idx); %For now, don't need to record TimeClosest
    Lat_Sample=Lat(idx);%Save latitude
    Lon_Sample=Lon(idx);%Save longitude
    Color_Sample=Nitr(j); %Change this to Nitr for nitrogen isotope
    
    %Make figure, plot track
    scatter(Lon_Sample,Lat_Sample,50,'black','o','filled');
    scatter(Lon_Sample,Lat_Sample,30,Color_Sample,'o','filled');
    
    end
    

end
 
%For main isotope whisker figure
figure; 
hold on
axis([0 50 Lon1 Lon2])

records = dir('*.mat');
for file = 1:max(size(records))
    
    load(records(file).name)     %Read in file

    %For each animal
    Lon=Track_Best(:,3);
    Lat=Track_Best(:,2);

    DownsampleNumber=max(round(size(Lon)/40))
    Lon_Sample=downsample(Lon,DownsampleNumber);
    Lat_Sample=downsample(Lat,DownsampleNumber);

    y=(0 + 4.*randn(max(size(Lat_Sample)),1));
    Color_Sample=Lon_Sample+y;

    %Make figure, plot track
    plot(1:max(size(Lon_Sample)),Lon_Sample,'LineWidth',0.1,'Color',[0.8 0.8 0.8])
    scatter(1:max(size(Lon_Sample)),Lon_Sample,50,'black','o','filled');
    scatter(1:max(size(Lon_Sample)),Lon_Sample,30,Color_Sample,'o','filled');
    
    clear Track_Best
end
 
 
 
 
%%%%%%%%%%%%%%%%% FOR ALL ANIMALS
Lon1=min(Lon)
Lon2=max(Lon)
Lat1=min(Lat)
Lat2=max(Lat)

