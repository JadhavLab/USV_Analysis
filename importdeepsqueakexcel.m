% Specify the folder where the files live.

myFolder = uigetdir;
while ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
end

% Get a list of all files in the folder with the desired file name pattern.

% you can use getAllFiles(dir,extension), this should help there with
% nested folders
filePattern = fullfile(myFolder, '*.xlsx'); % Change to whatever pattern you need.
theFiles = dir(filePattern);



%% Initialize table options
%{
opts = spreadsheetImportOptions("NumVariables", 16);

% Specify column names and types
opts.VariableNames = ["Call_ID", "Label", "Accepted", "Score", "BeginTimes", "EndTimes", "CallLengths", "PrincipalFrequencykHz", "LowFreqkHz", "HighFreqkHz", "DeltaFreqkHz", "FrequencyStandardDeviationkHz", "SlopekHzs", "Sinuosity", "MeanPowerdBHz", "Tonality"];
opts.VariableTypes = ["double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Accepted", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Accepted", "EmptyFieldRule", "auto");
opts.DataRange='A2';
%}

opts = spreadsheetImportOptions("NumVariables", 16);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:P1000";

% Specify column names and types
opts.VariableNames = ["ID", "Label", "Accepted", "Score", "BeginTimes", "EndTimes", "CallLengths", "PrincipalFrequencykHz", "LowFreqkHz", "HighFreqkHz", "DeltaFreqkHz", "FrequencyStandardDeviationkHz", "SlopekHzs", "Sinuosity", "MeanPowerdBHz", "Tonality"];
opts.VariableTypes = ["double", "categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Accepted", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Label", "Accepted"], "EmptyFieldRule", "auto");

%% Read in the first file to create a table to start with

combined_usvs = readtable(fullfile(theFiles(1).folder,theFiles(1).name), opts, "UseExcel", false);
combined_usvs(1,:)=[];
endtable=find(isnan(combined_usvs.ID),1,'first');
combined_usvs(endtable:end,:)=[];

%[pd, animal_id] = scrape_fileinfo(theFiles(1));

[cohort, animal_id,age] = scrape_fileinfo2(theFiles(1).name);
allRecordings=struct('rat',animal_id,'cohort',cohort,'age',age,...
    'datatable',combined_usvs,'sourcefile',fullfile(theFiles(1).folder,theFiles(1).name));

    



%% Loop through folder and add all other files to the big table

for k = 2:length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    [cohort, animal_id,age] = scrape_fileinfo2(theFiles(k).name);
    file_stats = readtable(fullFileName, opts, "UseExcel", false);
    file_stats(1,:)=[];
    endtable=find(isnan(file_stats.ID),1,'first');
    file_stats(endtable:end,:)=[];
    

    allRecordings(k).rat=animal_id;
    allRecordings(k).cohort=cohort;
    allRecordings(k).age=age;
    allRecordings(k).datatable=file_stats;
    allRecordings(k).sourcefile=fullFileName;
    
end

B = regexp(A,'\d*','Match')

%% lookup table for which animals are which genotype
% this is only for cohort 1
%{
lookuptable={'FF','FR','BR','BB','LL','BL','FL','RR';...
            'wt','het','fx','wt','wt','wt','het','wt';...
            'm', 'f', 'm', 'm', 'f', 'm', 'f', 'm';...
            '320','321','322','323','324','325','326','327'};
        
        % strcmpi strfind find contains 
for i=1:height(cohort2_full)
    tablematch=find(contains(lookuptable(1,:),cohort2_full.ratID(i,2:end)));
     cohort2_full.Genotype(i)=lookuptable(2,tablematch);
     cohort2_full.ratSex(i)=lookuptable(3,tablematch); 
     cohort2_full.ratNumber(i)=lookuptable(4,tablematch);  
     cohort2_full.cohort(i)=2;
end

        

%%

cohort1_lookuptable= {'1FL', '2FFLL', '1BL', '2FLBR', '2FRBL', '2LL', '1FR', '3FFLL', '2RR', '2BB', '1BR', '3FFRR';...
                      'fx', 'fx', 'wt', 'fx', 'wt', 'wt', 'fx', 'het', 'wt', 'wt', 'het', 'fx';...
                      'm',  'm', 'm', 'm', 'f', 'f', 'm', 'f', 'f', 'f', 'f', 'f';...
                      '301', '302', '303', '304','305', '306', '307','308', '309', '310', '311', '312'}; 
                  
 for i=1:height(cohort1_full)
     cellname= char(cohort1_full.ratID(i));
     tablematch=find(contains(cohort1_lookuptable(1,:), cellname));
     cohort1_full.Genotype(i)=cohort1_lookuptable(2,tablematch);
     cohort1_full.ratSex(i)=cohort1_lookuptable(3,tablematch);
     cohort1_full.ratNumber(i)=cohort1_lookuptable(4,tablematch);
     cohort1_full.cohort(i)=1; 
 end  
cohort1_full.Day=str2double(cohort1_full.Day);
                  
%% 
cohortfull=[cohort1_full; cohort2_full];
cohortfull(isnan(cohortfull.CallLengths),:)=[];
cohortfull.ratNumber=str2double(cohortfull.ratNumber);
%}

        
