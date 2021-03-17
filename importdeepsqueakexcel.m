% Specify the folder where the files live.

myFolder = uigetdir;
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.xlsx'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

%% Initialize table options
opts = spreadsheetImportOptions("NumVariables", 16);

% Specify column names and types
opts.VariableNames = ["Call_ID", "Label", "Accepted", "Score", "BeginTimes", "EndTimes", "CallLengths", "PrincipalFrequencykHz", "LowFreqkHz", "HighFreqkHz", "DeltaFreqkHz", "FrequencyStandardDeviationkHz", "SlopekHzs", "Sinuosity", "MeanPowerdBHz", "Tonality"];
opts.VariableTypes = ["double", "double", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Accepted", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Accepted", "EmptyFieldRule", "auto");

%% Read in the first file to create a table to start with

combined_usvs = readtable(theFiles(1).name, opts, "UseExcel", false);

[pd, animal_id] = scrape_fileinfo(theFiles(1));
    day_vec = repmat(pd,height(combined_usvs),1);
    id_vec = repmat(animal_id,height(combined_usvs),1);
    file_identifier_cols = array2table([day_vec, id_vec]);
combined_usvs = horzcat(file_identifier_cols, combined_usvs);
combined_usvs = combined_usvs(2:end,:);

%% Loop through folder and add all other files to the big table

for k = 2:length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    [pd,animal_id] = scrape_fileinfo(theFiles(k));

    file_stats = readtable(theFiles(k).name, opts, "UseExcel", false);
    day_vec = repmat(pd,height(file_stats),1);
    id_vec = repmat(animal_id,height(file_stats),1);
    
    file_identifier_cols = array2table([day_vec, id_vec]);

    file_stats = horzcat(file_identifier_cols, file_stats);
    combined_usvs = [combined_usvs; file_stats(2:end,:)];
end

combined_usvs