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

% now averaged across cohorts

% what do we want to plot
% lets toy with ncalls
allvars=fieldnames(allRecordings);
usevars=checkBox(allvars);
allvars=allvars(usevars);

cohortuse=1:3;
totext=@(a) cell2mat([a{1,1} {' '} a{1,2}]);

[rattypes,ua,ic]=unique(allRecs(:,{'sex','genotype'}),'rows');
for i=1:height(rattypes)
    rattypes.fullname{i}=cell2mat([rattypes{i,'sex'} ' ' rattypes{i,'genotype'}]);
end
ratuse=checkBox(rattypes.fullname);
    
mycolors=[1 0 .2; 1 .2 .0; 0 .2 1; .2 .0 1; 0.2 1 0; .0 1 .2];
for vr=1:length(allvars)

    rowuse=sum(ic==ratuse,2)>0;
    
    figure; 
    for k=1:length(ratuse)
    allX=allRecs.age(rowuse & ic==ratuse(k));
    allY=allRecs{rowuse & ic==ratuse(k),allvars{vr}}; 
    agemeans=accumarray(allX,allY,[],@nanmean);
    agevars=accumarray(allX,allY,[],@SEM);
    agemeans=agemeans(agemeans~=0);
    agevars=agevars(agemeans~=0);

    errorbar(2:2:length(agemeans)*2,agemeans,agevars,...
        'Color',mycolors(k,:),'LineWidth',2);
    hold on;
    end
    legend(rattypes.fullname(ratuse));
    ylabel(allvars{vr}); xlabel('Age');
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

[pcatable,ld,score]=pca(table2array(allRecs(:,9:end)));

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
        [~,uniquenaimals.peakln(i)]=max(animalsess.meanCallLengths);
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
    histogram(uniquenaimals.peakln(ic==i),1:.5:6.5);
    title(sprintf(' peak ln %s %s', ratGroup.sex{i}, ratGroup.genotype{i}));

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
[rattypes,ia,ic]=unique(allRecs(:,{'genotype','sex'}),'rows');

allRecs(ic==1,:)=[];

g1=table2array(allRecs(:,{'sex'}));
g2=table2array(allRecs(:,{'genotype'}));
g3=table2array(allRecs(:,{'age'}));

checked=checkBox(allRecs.Properties.VariableNames);
params=allRecs.Properties.VariableNames(checked);

for k=1:length(params)
    okdata=~isnan(allRecs.(params{k}));
    fprintf('\n %d %s \n',k,params{k});
    [a,b,c]=anovan(allRecs{okdata,params{k}}, {g1(okdata),g2(okdata),g3(okdata)},...
        'continuous',3,'model','full','sstype',2,...
        'varnames',{'sex','genotype','age'});
end

for k=1:length(params)
    okdata=~isnan(allRecs.(params{k})) & strcmpi(allRecs.sex,'m');
    fprintf('\n %d %s \n',k,params{k});
    [a,b,c]=anovan(allRecs{okdata,params{k}}, {g2(okdata),g3(okdata)},...
        'continuous',2,'model','full','sstype',3,...
        'varnames',{'genotype','age'});
    %multcompare(c);
    %keyboard
    %kill
end


[a,b,c]=anovan(allRecs.meanCallLengths,{g1,g2,g3},...
    'continuous',3,'model','full','sstype',2,...
    'varnames',{'sex','genotype','age'});

[a,b,c]=anovan(allRecs.meanFreq,{g1,g2,g3},...
    'continuous',3,'model','full','sstype',2,...
    'varnames',{'sex','genotype','age'});

[a,b,c]=anovan(allRecs.meanDeltaFreq,{g1,g2,g3},...
    'continuous',3,'model','full','sstype',2,...
    'varnames',{'sex','genotype','age'});

% observations


for k=1:length(params)
for i=2:height(rattypes)
    fprintf('rat type %s %s mean %s is %.3f +/- %.3f \n',...
        rattypes.genotype{i}, rattypes.sex{i},...
    params{k}, nanmean(allRecs{ic==i,params{k}}),...
    nanstd(allRecs{ic==i,params{k}}));
end
fprintf('\n');
end
        
        
%%
allX=allRecs.age;
allY=allRecs.meanCallLengths;
agemeans=accumarray(allX,allY,[],@nanmean);
agemeans=agemeans(agemeans~=0);
agevars=accumarray(allX,allY,[],@SEM);
agevars=agevars(agevars~=0);
figure; errorbar(2:2:length(agemeans)*2,agemeans,agevars);


%% now we find significant effects after accounting for age

% we'll report the difference between the two with effect size

for k=1:length(params)
    okdata=~isnan(allRecs.(params{k})) & strcmpi(allRecs.sex,'m');
    
    genodata=allRecs.genotype(okdata);
    age=allRecs.age(okdata);
    rawdata=allRecs{okdata,params{k}};
    
    agemeans=accumarray(age,rawdata,[],@mean);
    controlled=nan(length(rawdata),1);
    for i=1:length(agemeans)
        controlled(age==i)=rawdata(age==i)-agemeans(i);
    end
    
    % get the std by doign diff means/dprime means
    grpa=controlled(strcmpi(genodata,'fx'));
    grpb=controlled(strcmpi(genodata,'wt'));
    
    difference=nanmean(grpa)-nanmean(grpb);
    diffstd=difference/dprime(grpa,grpb);
    differencelow=difference-diffstd*2;
    
    fprintf('%s difference= %.4f +/- %.4f \n',params{k},difference,diffstd);
end


for k=1:length(params)
    okdata=~isnan(allRecs.(params{k})) & strcmpi(allRecs.sex,'m');
    
    genodata=allRecs.genotype(okdata);
    age=allRecs.age(okdata);
    rawdata=allRecs{okdata,params{k}};
    
    allages=unique(age);
    for i=1:length(allages)
        grpa=rawdata(age==allages(i) & strcmpi(genodata,'fx'));
        grpb=rawdata(age==allages(i) & strcmpi(genodata,'wt'));
        fprintf('%s age %d, diff=%.4f +/- %.4f \n',params{k},allages(i), nanmean(grpa)-nanmean(grpb),...
            2*(nanmean(grpa)-nanmean(grpb)/dprime(grpa,grpb)));
    end
end
