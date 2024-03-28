format compact
format long

%make "files" with TOPPIDs 
records = dir('*.mat');

%out=zeros(max(size(files)),3);

for file = 1:max(size(files))
    filename=records(file).name;
    t=load(filename);
    
    %out(file,1:3)=[t.TOPPID, t.MetaData.Group.CompleteTDR, t.MetaData.Group.CompleteTrack];
    if t.MetaData.Group.CompleteTDR==1
    
        out=horzcat(t.DiveLoc_Best(:,1:2),t.DiveStat,t.DiveType);
        writetable(out,[num2str(t.TOPPID) '_DiveStatType.csv'],'delimiter', ',');
        clear out
        
    else 
            
        out=(t.Track_Best);
        writetable(out,[num2str(t.TOPPID) '_Laura.csv'],'delimiter', ',');
        clear out
        
        
       
    end
    
    clear t

end

    