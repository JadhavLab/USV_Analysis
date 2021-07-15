% full plotanalysis

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

param{1}='PCA2';
param{2}='PCA1';

% first only pull males

allmales=allRecs(strcmpi(allRecs.sex,'m'),:);
    [days,unidx,idx]=unique(allmales.age);
    mycolors=parula(length(unidx)+2);
    figure;
for i=1:height(allmales)
    %if days(idx(i))<10
    if strcmpi(allmales.genotype(i),'wt')
        hs(1)=scatter(allmales{i,param{1}},allmales{i,param{2}},14,mycolors(idx(i),:)); hold on;
    else
        hs(2)=scatter(allmales{i,param{1}},allmales{i,param{2}},14,mycolors(idx(i),:),'filled'); hold on;
    end
    %end
end
    legend(hs,{'WT','FX'});    

    
% lets run some ancovas

[pcatable,ld,score]=pca(table2array(allRecs(:,9:16)));

allRecs.PCA1=ld(:,1);
allRecs.PCA2=ld(:,2);
allRecs.PCA3=ld(:,3);
%% now some questions
%{
1/ do se see any group effects? maybe run an ancova on *any* output across
genotype day and sex

% what if you pca these animals, do you see any groupings?

% you could pca the calls themselves somehow to cluster
%}

% cut the data in half
PCAsort=sort(allRecs.PCA1(~isnan(allRecs.PCA1)));
PCAind=elbow_method(PCAsort,[],[],false);
PCAthresh=PCAsort(PCAind);

% now peek which day each animal achieves this
[uniqueanimals,ia,ic]=unique(allRecs(:,[2 6]));

for i=1:height(uniqueanimals)
    animalsess=sortrows(allRecs(ic==i,:),3);
    try
        uniqueanimals.daycrit(i)=find(animalsess.PCA1<PCAthresh,1,'last');
    catch
        uniqueanimals.daycrit(i)=nan;
    end
    uniqueanimals.sex(i)=allRecs.sex(ia(i));
    uniqueanimals.genotype(i)=allRecs.genotype(ia(i));
end

[ratGroup,ia,ic]=unique(uniqueanimals(:,[4 5]));
figure;
for i=1:height(ratGroup)
    ratGroup.daysCrit(i)={uniqueanimals.daycrit(ic==i)};
    subplot(height(ratGroup),1,i);
    histogram(uniqueanimals.daycrit(ic==i),1:2:17);
    title(sprintf('%s %s', ratGroup.sex{i}, ratGroup.genotype{i}));
end


%%


figure;
for i=1:height(ratGroup)
    ratGroup.daysCrit(i)={uniqueanimals.daycrit(ic==i)};
    subplot(height(ratGroup),1,i);
    histogram(uniqueanimals.daycrit(ic==i),1:2:17);
    title(sprintf('%s %s', ratGroup.sex{i}, ratGroup.genotype{i}));
end

%% anovans
g1=table2array(allRecs(:,{'sex'}));
g2=table2array(allRecs(:,{'genotype'}));
g3=table2array(allRecs(:,{'age'}));
[a,b,c]=anovan(allRecs.meanCallLengths,{g1,g2,g3},...
    'continuous',3,'model','full','sstype',1,...
    'varnames',{'sex','genotype','age'});

[a,b,c]=anovan(allRecs.meanFreq,{g1,g2,g3},...
    'continuous',3,'model','full','sstype',1,...
    'varnames',{'sex','genotype','age'});

[a,b,c]=anovan(allRecs.meanTonality,{g1,g2,g3},...
    'continuous',3,'model','full','sstype',1,...
    'varnames',{'sex','genotype','age'});
