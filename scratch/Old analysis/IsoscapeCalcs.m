%%%%%%%%%%%%%%
%CALCULATE METRICS FROM MAT FILES

filenames=dir('*.mat');
out=zeros(max(size(filenames)),3);

for loop=1:size(filenames,1)
    
    load(filenames(loop).name)
   


    %TOPPID
    out(loop,1)=TOPPID;
    
    %maxlat
    out(loop,2)=ForagingSuccess.TBGEGainPercent;

    
end

%%%%%%%%%%%%%%
%CALCULATE METRICS FROM DRIFT ANALYSIS MAT FILES

filenames=dir('*.mat');
out=zeros(max(size(filenames)),4);

for loop=1:size(filenames,1)

    load(filenames(loop).name)
    
    out(loop,1)=TOPPID;
    
    %MeanDriftRate
    out(loop,2)=mean(DriftDives4(:,2));
    
    %SumDriftRate
    out(loop,3)=sum(DriftDives4(:,2));
    
    %Date of switch to positive buoyancy
    out(loop,4)=DriftDives4(find(DriftDives4(:,2)==max(DriftDives4(:,2))),1);
    
end