% This code is for importing the tables Track_Best (has lat and long), 
% DiveType (has drift rates for each dive), and DriveStat (has dates for each dive)
% from the MATLAB files of several seals and export them as CSV files

% Written by: FV, JJ
% Project: Eseals Isoscapes
% Last modified: 31 Jan 2024

% Retrieve from the directory all the dive MATLAB files whose TOPPID 
% correspond to years later than 2010 and put them into a list:
listing=dir("202*.mat");

for x=1:size(listing,1) %the itineration goes from the first MATLAB file
    %until the length of the "listing" table 

    % Loading column "names" that exists in the table "listing", which is
    % the name of the MATLAB file:
    load(listing(x).name);
    
    if size(Track_Best,1)~=0 %If Track_Best has a number of rows that is different from
        %zero, go on with the for loop; otherwise, do nothing (this
        %prevents the code from stopping if for any given seal, there's no
        %Track_Best table)
        
        % Add a column for Seal ID to TrackBest
        Track_Best.SealID=repmat(MetaData.FieldID,size(Track_Best,1), 1); 
        % Creating a matrix that will have a size according to the number 
        % of rows in table Track_Best, and one column. Function "size" determines 
        % the size of the matrix. The matrix is populated with the FieldID 
        % (which is the Seal ID) within the table MetaData. 
        
        % Export the table Track_Best into a CSV file
        writetable(Track_Best(:,{'SealID', 'JulDate', 'Lat', 'Lon'}), ...
            strjoin(["Track_Best_",TOPPID,'_',MetaData.FieldID,'.csv'],""));
            % To extract the table Track_Best as a CSV files, with the specified 
            % columns, and then specifying the CSV file name
    end %this ends the "if" statement

    if size(DiveType,1)~=0 %If DiveType has a number of rows that is different from
        %zero, go on with the for loop; otherwise, do nothing (this
        %prevents the code from stopping if for any given seal, there's no
        %DiveType table)
        
        % Add a column for Seal ID to DiveType
        DiveType.SealID=repmat(MetaData.FieldID,size(DiveType,1), 1); 
        % Creating a matrix that will have a size according to the number 
        % of rows in table DiveType, and one column. Function "size" determines 
        % the size of the matrix. The matrix is populated with the FieldID 
        % (which is the Seal ID) within the table MetaData. 
        
        % Export the table DiveType into a CSV file
        writetable(DiveType(DiveType.DiveType==2,{'SealID', 'DiveNumber', ...
            'DriftRate', 'DiveType', 'DiveTypeName'}), ...
            strjoin(["DiveType_",TOPPID,'_',MetaData.FieldID,'.csv'],""));
            % To extract the table DiveType as a CSV files, containing only the 
            % rows with DiveType values of 2 and the named columns, and then 
            % specifying the CSV file name

    end %this ends the "if" statement

     if size(DiveStat,1)~=0 %If DriveStat has a number of rows that is different from
        %zero, go on with the for loop; otherwise, do nothing (this
        %prevents the code from stopping if for any given seal, there's no
        %DiveStat table)
        
        % Add a column for Seal ID to DriveStat
        DiveStat.SealID=repmat(MetaData.FieldID,size(DiveStat,1), 1); 
        % Creating a matrix that will have a size according to the number 
        % of rows in table DriveStat, and one column. Function "size" determines 
        % the size of the matrix. The matrix is populated with the FieldID 
        % (which is the Seal ID) within the table MetaData. 
        
        % Export the table DriveStat into a CSV file
        writetable(DiveStat(:,{'SealID', 'DiveNumber', 'Year', 'Month','Day', ...
            'Hour','Min','JulDate'}), ...
            strjoin(["DiveStat_",TOPPID,'_',MetaData.FieldID,'.csv'],""));
            % To extract the table DiveStat as a CSV files, with the specified 
            % columns, and then specifying the CSV file name

    end %this ends the "if" statement

end

 

