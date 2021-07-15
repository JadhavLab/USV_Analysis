% preprocess the data
%
%
%

clearvars -except allRecordings
load('USVmetadata.mat');

%% now get some summary stats for each rat
USVcohorts.LegMarkings=strtrim(string(USVcohorts.LegMarkings));
USVcohorts.ToeMarkings=strtrim(string(USVcohorts.ToeMarkings));


for i=1:length(allRecordings)
    
    % first line up the animal characteristics
    mymatch= find(str2double(allRecordings(i).cohort(2))==USVcohorts.Cohort & ...
    (strcmpi(allRecordings(i).rat,USVcohorts.Animal) | ...
    strcmpi(allRecordings(i).rat,USVcohorts.LegMarkings) |...
    strcmpi(allRecordings(i).rat,USVcohorts.ToeMarkings)));
    if ~isempty(mymatch)

    allRecordings(i).animal=char(USVcohorts.Animal(mymatch));
    allRecordings(i).genotype=char(USVcohorts.Genotype(mymatch));
    allRecordings(i).sex=char(USVcohorts.Sex(mymatch));
    allRecordings(i).cohort=str2double(allRecordings(i).cohort(2));
    allRecordings(i).age=str2double(allRecordings(i).age(2:end));
    end
end


for i=1:length(allRecordings)
    % now summary statistics
    allRecordings(i).meanCallLengths=nanmean(allRecordings(i).datatable.CallLengths);
    allRecordings(i).varCallLengths=nanvar(allRecordings(i).datatable.CallLengths);
    allRecordings(i).meanFreq=nanmean(allRecordings(i).datatable.PrincipalFrequencykHz);
    allRecordings(i).varFreq=nanvar(allRecordings(i).datatable.PrincipalFrequencykHz);
    allRecordings(i).meanDeltaFreq=nanmean(allRecordings(i).datatable.DeltaFreqkHz);
    allRecordings(i).varDeltaFreq=nanvar(allRecordings(i).datatable.DeltaFreqkHz);

    allRecordings(i).meanFreqSTD=nanmean(allRecordings(i).datatable.FrequencyStandardDeviationkHz);
    allRecordings(i).meanPowerdBHz=nanmean(allRecordings(i).datatable.MeanPowerdBHz);
    allRecordings(i).meanTonality=nanmean(allRecordings(i).datatable.Tonality);
    end


allRecs=struct2table(allRecordings);
    
%%
load Cohorts1-4.mat

