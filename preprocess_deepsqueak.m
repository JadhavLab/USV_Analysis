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

%% now view


% what do we want to plot
% lets toy with ncalls
allvars=fieldnames(allRecordings);
usevars=checkBox(allvars);
allvars=allvars(usevars);

for vr=1:length(allvars)

thisvar=allvars{vr};

% for each cohort
figure;
for ch=1:4
    thisCohort=allRecs(allRecs.cohort==ch,:);
    %thisCohort.animal=cell2mat(thisCohort.animal);
    % now plot each bit by age...
    % we want something like rat by age
    allrats=unique(thisCohort.animal);
    alldays=unique(thisCohort.age);
    
    varmat=[]; idxmat=[];
    [genomatch,unidx,idx]=unique(thisCohort(:,[7 8]),'rows');
    colors=lines(length(unidx));
    subplot(1,4,ch);
    for i=1:length(allrats)
        for j=1:length(alldays)
            matmatch=find(strcmpi(thisCohort.animal,allrats(i))...
                & thisCohort.age==alldays(j));
            if ~isempty(matmatch)
            varmat(i,j)=thisCohort.(thisvar)(matmatch);
            idxmat(i,j)=idx(matmatch);
            else
                varmat(i,j)=nan; idxmat(i,j)=nan;
            end
        end
        plot(alldays,varmat(i,:),'Color',colors(idxmat(i,j),:));
        hold on;
    end
    title(sprintf('%s, cohort %d',thisvar,ch));
end
end
%% and not by cohort


for vr=1:length(allvars)

thisvar=allvars{vr};

% for each cohort
figure;

    %thisCohort.animal=cell2mat(thisCohort.animal);
    % now plot each bit by age...
    % we want something like rat by age
    allrats=unique(allRecs.animal);
    alldays=unique(allRecs.age);
    
    varmat=[]; idxmat=[];
    [genomatch,unidx,idx]=unique(allRecs(:,[7 8]),'rows');
    for i=1:height(genomatch)
        unnames{i}=[cell2mat(genomatch{i,1}) ' ' cell2mat(genomatch{i,2})];
    end
    mycolors=(jet(length(unidx)+2));

    for i=1:length(allrats)
        for j=1:length(alldays)
            matmatch=find(strcmpi(allRecs.animal,allrats(i))...
                & allRecs.age==alldays(j));
            if ~isempty(matmatch)
            	varmat(i,j)=allRecs.(thisvar)(matmatch);
            	idxmat(i,j)=idx(matmatch);
            else
                varmat(i,j)=nan; idxmat(i,j)=nan;
            end
        end
        hp(i)=plot(alldays,varmat(i,:),'Color',mycolors(mode(idxmat(i,:)),:));
        hold on;
    end
    legend(hp(length(hp)-unidx),unnames);
    title(sprintf('%s All Cohorts',thisvar));
    
end
%% lets plot two, color them, open fx males, filled control males

param{1}='meanFreq';
param{2}='meanTonality';

% first only pull males

allmales=allRecs(strcmpi(allRecs.sex,'m'),:);
    [days,unidx,idx]=unique(allmales.age);
    mycolors=parula(length(unidx)+2);
    figure;
for i=1:height(allmales)
    if idx(i)<7
    if strcmpi(allmales.genotype(i),'wt')
        hs(1)=scatter(allmales{i,param{1}},allmales{i,param{2}},14,mycolors(idx(i),:)); hold on;
    else
        hs(2)=scatter(allmales{i,param{1}},allmales{i,param{2}},14,mycolors(idx(i),:),'filled'); hold on;
    end
    end
end
    legend(hs,{'WT','FX'});    

        



%% now some questions
%{
1/ do se see any group effects? maybe run an ancova on *any* output across
genotype day and sex

% what if you pca these animals, do you see any groupings?

% you could pca the calls themselves somehow to cluster
%}