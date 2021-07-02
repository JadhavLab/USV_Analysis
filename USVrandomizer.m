%% import the data

opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Summary";
opts.DataRange = "A2:E45";

% Specify column names and types
opts.VariableNames = ["Cohort", "Animal", "Genotype", "Sex", "markings"];
opts.VariableTypes = ["double", "string", "categorical", "string", "string"];

% Specify variable properties
opts = setvaropts(opts, ["Animal", "Sex", "markings"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Animal", "Genotype", "Sex", "markings"], "EmptyFieldRule", "auto");

% Import the data
USVrecordings = readtable("/Users/johnbladon/Downloads/USV recordings .xlsx", opts, "UseExcel", false);

% clear variables
clear opts
%% create the age matrix

USVrecordings(isnan(USVrecordings.cohort),:)=[];

% now blow it up for early mid and late

%% now get the unique categories

%% now pull even numbers from each

%% we need ~4*8 sessions


%% pull all the files

mydir=uigetdir();

allfiles=getAllFiles(mydir,'.wav');
