

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "Summary";
opts.DataRange = "A1:G45";

% Specify column names and types
opts.VariableNames = ["Cohort", "Animal", "Genotype", "Sex", "EarMarkings", "ToeMarmings", "LegMarkings"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "string", "string"];

% Specify variable properties
opts = setvaropts(opts, ["Animal", "EarMarkings"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Animal", "Genotype", "Sex", "EarMarkings", "ToeMarmings", "LegMarkings"], "EmptyFieldRule", "auto");

% Import the data
USVcohorts = readtable("C:\Users\Jadhavlab\Desktop\USV recordings .xlsx", opts, "UseExcel", false);


% Clear temporary variables
clear opts

% clean table up a bit

USVcohorts(isnan(USVcohorts.Cohort),:)=[];
USVcohorts.Sex=lower(string(USVcohorts.Sex));
USVcohorts.Properties.VariableNames{6} = 'ToeMarkings';


%% now grab all files

mydir=uigetdir;

subdirs=dir(mydir); subdirs(1:2)=[];

USVfiles=dir(fullfile(subdirs(1).folder,subdirs(1).name));
USVfiles(1:2)=[];
for k=2:length(subdirs)
    tmpfiles=dir(fullfile(subdirs(k).folder,subdirs(k).name));
    tmpfiles(1:2)=[];
    USVfiles=[USVfiles; tmpfiles];
end
clear tmpfiles;
    

%% now go through each file and assign a cohort, a day, and an animal


for i=1:length(USVfiles)
    
    nameseg=find(USVfiles(i).name=='_');
    
    if length(nameseg)<2
        nameseg=find(USVfiles(i).name=='-');
    end
    
    USVfiles(i).Cohort=USVfiles(i).name(1:nameseg(1)-1);
    USVfiles(i).Day=USVfiles(i).name(nameseg(1)+1:nameseg(2)-1);
    USVfiles(i).rat=USVfiles(i).name(nameseg(2)+1:end-4);
end

%% now match up the animals with their xf name, sex and genotype


for i=1:length(USVfiles)
    thisCohort=USVcohorts(USVcohorts.Cohort==str2double(USVfiles(i).Cohort(2:end)),:);
    thisRat=thisCohort(strcmpi(thisCohort.ToeMarkings,USVfiles(i).rat),:);
    if height(thisRat)==0
        thisRat=thisCohort(strcmpi(thisCohort.LegMarkings,USVfiles(i).rat),:);
    end
    if height(thisRat)==0
        thisRat=thisCohort(strcmpi(thisCohort.Animal,USVfiles(i).rat),:);
    end
    if height(thisRat)>0
        USVfiles(i).Rat=thisRat.Animal;
        USVfiles(i).Genotype=thisRat.Genotype;
        USVfiles(i).Sex=thisRat.Sex;
        USVfiles(i).EarMarkings=thisRat.EarMarkings;
        USVfiles(i).ToeMarkings=thisRat.ToeMarkings;
        USVfiles(i).LegMarkings=thisRat.LegMarkings;
    else
        fprintf('USV file %s had no corresponding rat \n', USVfiles(i).name);
    end
end

    
    
    
    