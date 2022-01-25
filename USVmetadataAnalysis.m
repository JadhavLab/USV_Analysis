%% these data come from a webpage that you must download, heres the page

web('https://docs.google.com/spreadsheets/d/1NQK8Sv2n_1EUyIOZScPLmUXlLDKu54r3QYa5IO0BRcE/edit#gid=1940592088');

% then download the whole doc, and then use that to pull the below
%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
opts.Sheet = "Summary";
opts.DataRange = "A1:H200";

% Specify column names and types
opts.VariableNames = ["cohort", "animal", "genotype", "sex", "earMarkings", "toeMarkings", "legMarkings","DOB"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "string", "string", "datetime"];

% Specify variable properties
opts = setvaropts(opts, ["animal", "earMarkings"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["animal", "genotype", "sex", "earMarkings", "toeMarkings", "legMarkings","DOB"], "EmptyFieldRule", "auto");

% Import the data
try
USVcohorts = readtable("C:\Users\Jadhavlab\Desktop\USV recordings.xlsx", opts, "UseExcel", false);
catch
    USVcohorts = readtable("C:\Users\John Bladon\Desktop\USV recordings.xlsx", opts, "UseExcel", false);
end

% Clear temporary variables
clear opts

% clean table up a bit
USVcohorts(isnan(USVcohorts.cohort),:)=[];
USVcohorts.sex=lower(string(USVcohorts.sex));



% now import the metadata for each recording
% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
Sheetnames=["Cohort 1","Cohort 2","Cohort 3","Cohort 4","Cohorts 5 & 6","Cohort 7","Cohort 8","Cohort 9"];

for i=1:length(Sheetnames)
    opts.Sheet = Sheetnames(i);
    opts.DataRange = "A3:H300";
    
    % Specify column names and types
    opts.VariableNames = ["Date", "FileName", "Animal", "Weight", "Time", "XFName", "OtherID", "VarName8"];
    opts.VariableTypes = ["datetime", "string", "string", "double", "string", "string", "string", "string"];
    
    % Specify variable properties
    opts = setvaropts(opts, ["FileName", "Animal", "Time", "XFName", "OtherID", "VarName8"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["FileName", "Animal", "Time", "XFName", "OtherID", "VarName8"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "Date", "InputFormat", "");
    
    % Import the data
    CohortsTemp = readtable("C:\Users\John Bladon\Desktop\USV recordings.xlsx", opts, "UseExcel", false);
    CohortsTemp(CohortsTemp.FileName=="",:)=[];
    if i==1
        USVSession=CohortsTemp;
    else
        USVSession=vertcat(USVSession,CohortsTemp);
    end
end

USVSession(cellfun(@(a) isempty(a),USVSession.Animal),:)=[];
% Clear temporary variables
clear opts

% now clean up
for i=1:height(USVSession)
    oldseg=find(USVSession.FileName{i}=='-');
    USVSession.FileName{i}(oldseg)='_';
    
    nameseg=find(USVSession.FileName{i}=='_');
    USVSession.cohort(i)=str2double(USVSession.FileName{i}(2:nameseg(1)-1));
    USVSession.day(i)=str2double(USVSession.FileName{i}(nameseg(1)+2:nameseg(2)-1));
    % and we hold off on animal name, we already have quite a few names
    % here
end
%%
%
%     wavfile time!!!!
%
%
%
%
%
%
% now grab all WAVfiles

% i keep them in a single folder
[wavfile,wavdir]=uigetfile('make sure you have all the wav files in this folder');

USVfiles=dir(wavdir);
USVfiles(1:2)=[];

    

% now go through each file and assign a cohort, a day, and an animal


for i=1:length(USVfiles)
    % first convert all dashes to underscores to standardize
    oldseg=find(USVfiles(i).name=='-');
    USVfiles(i).name(oldseg)='_';

    nameseg=find(USVfiles(i).name=='_');

    
    USVfiles(i).cohort=str2double(USVfiles(i).name(2:nameseg(1)-1));
    USVfiles(i).day=str2double(USVfiles(i).name(nameseg(1)+2:nameseg(2)-1));
    USVfiles(i).rat=USVfiles(i).name(nameseg(2)+1:end-4);
end

% sort by recording date
[~,index] = sortrows([USVfiles.datenum].'); USVfiles = USVfiles(index); clear index

USVfiles=struct2table(USVfiles);
%% now match up the animals with their xf name, sex and genotype

% It will be USV recordings.Session, it will have:
% 1. date
% 2. Filename
% 3. directory
% 4. Animal (XF)
% 5. Weight
% 6. Age
% 7. Sex
% 8. Genotype
% 9. Cohort
% 10. Other markings/aliases
% we'll tack the cohort info to each wavfile

% USVSession has the most entries, so we'll use that, and slot in variables
% there

    allSessions=table({},{},{},[],[],[],{},[],[],'VariableNames',...
    {'fileName','folder','date','datenum','cohort','age','aliases','weight','XFname'});
for i=1:height(USVSession)
     % first grab the data from the filename
    easymatch=USVfiles.cohort==USVSession.cohort(i) & USVfiles.day==USVSession.day(i);

    mygroup=USVfiles(easymatch,:);
    
    % now try matching animal- we will match rat in mygroup to any in the
    % XFName, OtherID, or Animal
    matches=[strcmpi(USVSession.Animal(i),mygroup.rat) strcmpi(USVSession.OtherID(i),mygroup.rat)];

    [nmatches,matchinds]=sort(sum(matches,2), 'descend');

    if ~isempty(matches) && (nmatches(1)>0)
        fprintf('  %s matches %s or %s %d times, using this \n', mygroup.rat{matchinds(1)},...
            USVSession.Animal(i),USVSession.OtherID(i), nmatches(1))
        % fill in from the filedata first
        allSessions(i,[1:7])=mygroup(matchinds(1),[1 2 3 6 7 8 9]);

    else
        fprintf('...found no matches for %s, moving on...\n',USVSession.Animal(i));
        % else fill in from your metadata
        allSessions.fileName{i}=USVSession.FileName{i};
        allSessions.date{i}=datestr(USVSession.Date(i));
        allSessions.cohort(i)=USVSession.cohort(i);
        allSessions.age(i)=USVSession.day(i);
    end
    % these always come from metadata
    allSessions.aliases{i,:}=[USVSession.Animal(i) USVSession.OtherID(i)];
    allSessions.weight(i)=USVSession.Weight(i);
    allSessions.XFname(i)=str2double(USVSession.XFName{i}(3:end));
    if i==1
        w = warning('query','last');
        id = w.identifier;
        warning('off',id)
    end
end
    warning('on',id)

%% this is just chaff really

[a,b,c]=unique(USVfiles(:,[7 11 12 16]),'rows');


Analysts={'Norelis','Jay','Zach'};

for i=1:height(a)
    % cycle through analysts (this ougth to be pretty uniform)
    a.Analyst{i}=Analysts{mod(i,3)+1};
    % gather the indices for our source pool
    sourcepool=find(c==i);
    % now pick a random session from that pool
    sourceind=sourcepool(randi(length(sourcepool),1));
    a.sourcefile{i}=USVfiles.name(sourceind);
    a.fileloc{i}=USVfiles.folder(sourceind);
end

USVpreprocessLedger=a;
% so we have 32 unique categories, so now we have to distribute them across
% the experimenters

% now divvy:
%{
So i want to balance across sex, genotype, cohort, and early v late.
  
The mid point of recordings is P10.  We can justify this because after P20,
calls fall off in frequency
  
%}
    
    