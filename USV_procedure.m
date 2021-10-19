%{
This is the analysis sketchpad of the usv data from start to finish
 I will start with blocks of sketch code here and then offload them into
 functions or other scripts if they become useful.

the preprocessing pipeline was as follows:
-data were acquired, see acquisition protocol
-a systematically random selection of files were pulled and manually
curated for training data.
-a network was built on those calls using DeepSqueak v 3.0.1
-the full detection algorithm was then performed on the full gambit of data
using both the rat detector yolo r1 and rat detector JHB networks together
-a post-hoc denoiser was then applied to ALL FILES

Manual edits:
First some observations:

-the bounding boxes tend to miss the start and end of the calls, leaving
loud bits outside the box

- The contour isnt quite right sometimes, need to check how they threshold
to get a contour.
    The contour math:
    -they use a conjunctive approach.  they need at least 6 x coords in
    which... both the entropy (geomean/mean)>a threshold) and the loudness
    threshold are used. Loudness threshold=.825 and entropy is .215.. then
    they increase each by 10% until you get 6 dots
    -the elbow method could be used to find the top bit of spiking
    -you could also use a slope, e.g. when does the area between the linear
    regression and the real data increase by less than say 10%

-effectively, I think we can use a 'loudness' threshold because thats
effectively what a 'listening mother' will actually hear, and therefore is
whats functionally relevant
- we can ignore the contour, it appears that uvseg has a very easy and
intuitive way of grabbing contours.  SO the to do list:

1. reject boxes that dont have enough 'intensity'- this will be a batch
process that modfies the output files that deepsqueak makes
    a. to do this first get a spectrogram of the full day, i think there
    is an option for median, and then use that as a background subtract for
    all pixels.
    b. now get an average peak intensity for all calls in that day, and
    subtract anything less than say 3sd below that, or like less than say
    50% of the distance from that to the mean.
    c. assign all those calls 'reject' calls, or delete those rows

2. regather the contours based on the usvseg method, and then see whether
we need to combine boxes or expand boxes.

3. use uvseg to gather the principal contour(s) and analyze those for
frequency- maybe look at principal frequencie(s) using a findpeaks method,
and then run regressions on those lines to get pitch shifts vs slides

4. use a machine learning algorithm, probably using the full column to
segment usvs, dont fuck with it, it should just work.



%}

%% run this procedure to gather ROIs and remove squeaks that are too soft

% one remaining problem is that I still have low frequency noise being
% detected as signal.  I think i may be able to get over this by running it
% over a gradient smoothing algorithm.  the idea is that there should be a
% high level of at least somewhat locally high gradients that the noise
% wont have.  I can probably get around this by measuring the skew (geomean
% over mean) over say the 10,000 hz surrounding the spline.  the spline
% also shouldnt be too scattered...

edit removeSoftCalls

try
    load('USVmetadata2021-10-15.mat');
catch
    edit USVmetaAggregator 
end

edit mergeMetaCalls

%% just curious about weight
weightmat={}; allweights={};
for i=1:max(USVSession.Cohort)
    uniquedays=unique(USVSession.Age(USVSession.Cohort==i));
    allweights{i}=uniquedays;
    for j=1:length(uniquedays)
        weightmat{i}(:,j)=USVSession.Weight(USVSession.Cohort==i & USVSession.Age==uniquedays(j));
    end
end
figure;
clear eb;
for i=1:length(weightmat)
    eb(i)=errorbar(allweights{i},nanmean(weightmat{i}),nanstd(weightmat{i}));
    hold on;
end
xlabel('Age pnd'); ylabel('Weight, gm')
legend(eb)
% it appears that cohort 5 is older than reported...
%%

load('C:\Users\John Bladon\Desktop\USVdataFULL2021-10-15.mat');
% this is what deepsqueak gives me to work with...

%% now that we have a datatable, lets see what falls out
% first lets plot everything by sex/genotype and time
oktouse=find(cellfun(@(a) ~isempty(a), USVSession.CallStats(:)));
useSession=USVSession(oktouse,:);
allstats=useSession.CallStats{oktouse(1)}.Properties.VariableNames;
%usestats=allstats(checkBox(allstats,'Choose the stats you want to see'));
usestats=allstats(7:17); %just use all these ones...
% for each usestat
% for each session
% gather the average stat
% errorbar for each genotype
% back color plot for each cohort of that genotype

[genotypes,~,idx]=unique(useSession(:,{'Sex','Genotype'}),'rows');
genonames={'het F','wt F','fx M','wt M'};
gtcolors=lines(4);

mystat=cellfun(@(a) height(a),useSession.CallStats(:));
figure; 
sp=subplot(3,4,1);
clear eb
% this is call counts
for gt=1:4 % 4 genotypes
    daymean=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@mean,nan);
    daydev=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@std,nan);
    alldays=unique(useSession.Age(idx==gt));
    eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
        'Color',gtcolors(gt,:));
    hold on;
end
xlabel('Postnatal Age');
title(sprintf('Num calls'));
legend(eb,genonames)
% this is all the other parameters
for i=1:length(usestats)
    mystat=cellfun(@(a) mean(a.(usestats{i}),'omitnan'),useSession.CallStats(:));
    sp(i+1)=subplot(3,4,i+1);
    for gt=1:4 % 4 genotypes
        daymean=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@mean,nan);
        daydev=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@std,nan);
        alldays=unique(useSession.Age(idx==gt));
        eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
            'Color',gtcolors(gt,:));
        hold on;
    end
    xlabel('Postnatal Age'); ylabel(sprintf('Mean %s',usestats{i}));
    title(sprintf('Mean %s',usestats{i}));
    legend(eb,genonames)
end

linkaxes(sp,'x'); xlim([3 21]);

% now for the call type analysis, which I think will pull something out
% here
%% now plot these against weight
% are fx animals heavier than their counterparts?

figure; sp=subplot(1,3,1);

for gt=1:4 % 4 genotypes
    daymean=accumarray(useSession.Age(idx==gt),useSession.Weight(idx==gt),[],@nanmean,nan);
    daydev=accumarray(useSession.Age(idx==gt),useSession.Weight(idx==gt),[],@nanstd,nan);
    alldays=unique(useSession.Age(idx==gt));
    eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
        'Color',gtcolors(gt,:));
    hold on;
end
linkaxes(sp,'x'); xlim([3 21]);
ylabel('Weight, gms'); xlabel('Age, pnd');

% lets see if this is real
[p,tbl,stats]=anovan(useSession.Weight,{useSession.Age, idx},'continuous',1,...
    'varnames',{'Age','Genotype'},'model','full');

[p,tbl,stats]=anovan(useSession.Weight,{useSession.Age, idx==3},'continuous',1,...
    'varnames',{'Age','Genotype'});


%% so I think we should plot the p(call) for each call type for each animal

allcalls=cellfun(@(a) a.Label, useSession.CallStats(:), 'UniformOutput', false);
everycall=cell2mat(cellfun(@(a) str2double(a), allcalls, 'UniformOutput', false));
everycall(isnan(everycall))=26; % there are 26 types of call (26 is unidentified)
uniquecalls=unique(everycall);
for i=1:length(uniquecalls)
    figure;
    mystat=cellfun(@(a) mean(str2double(a.Label)==i,'omitnan'),useSession.CallStats(:));
    %sp=subplot(3,4,i+1);
    for gt=1:4 % 4 genotypes
        daymean=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@nanmean,nan);
        daydev=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@nanstd,nan);
        alldays=find(~isnan(daymean));
        eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
            'Color',gtcolors(gt,:));
        hold on;
    end
    % pull out random images for those clusters
    %sp(2)=
    xlabel('Postnatal Age');
    title(sprintf('call %d',i));
    legend(eb,genonames)
end

%% howabout call bouts

% where is the peak in the i-c-i histogram?
% two ways of doing this: create a sparse mat and xcorr it
% or histogram the ISIs

% maybe the peak of that isi histogram?
xedges=0:.005:.4;
xcenters=mean([xedges(1:end-1); xedges(2:end)]);
isiPeak=[];
callDurPeak=[];
wb=waitbar(0,'starting');
for i=1:height(useSession)
    allISI=diff([useSession.CallStats{i}.EndTimes(1:end-1) useSession.CallStats{i}.BeginTimes(2:end)],1,2);
    [y,x]=histcounts(allISI,0:.005:.4);
    [~,peakind]=max(y); isiPeak(i)=xcenters(peakind);
    [y,x]=histcounts(useSession.CallStats{i}.CallLengths,xedges);
    [~,peakind]=max(y); callDurPeak(i)=xcenters(peakind);
    waitbar(i/height(useSession),wb,'running now');
end
close(wb)
figure;
%sp=subplot(3,4,i+1);
for gt=1:4 % 4 genotypes
    daymean=accumarray(useSession.Age(idx==gt),isiPeak(idx==gt),[],@nanmean,nan);
    daydev=accumarray(useSession.Age(idx==gt),isiPeak(idx==gt),[],@nanstd,nan);
    alldays=find(~isnan(daymean));
    eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
        'Color',gtcolors(gt,:));
    hold on;
end
xlabel('Postnatal Age');
ylabel('Peak inter-call-interval, seconds')
legend(eb,genonames)

% how many calls are within say 100 msec of eachother?
figure;
%sp=subplot(3,4,i+1);
for gt=1:4 % 4 genotypes
    daymean=accumarray(useSession.Age(idx==gt),callDurPeak(idx==gt),[],@nanmean,nan);
    daydev=accumarray(useSession.Age(idx==gt),callDurPeak(idx==gt),[],@nanstd,nan);
    alldays=find(~isnan(daymean));
    eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
        'Color',gtcolors(gt,:));
    hold on;
end
xlabel('Postnatal Age');
ylabel('Call Length Mode, seconds')
legend(eb,genonames)


% ncalls within n milliseconds

% 

%% cluster early day animals, cluster late day animals





%% plot certain characteristics against eachother

[alldays,~,dayind]=unique(useSession.Age);
daycolor=parula(length(alldays));
daycolor=daycolor(dayind,:);

genocolor=gtcolors(idx,:);

Param1='CallLengths';
Param2='Sinuosity';


usestats=allstats(7:17); %just use all these ones...
for i=1:length(usestats)-1
    for j=i:length(usestats)
        Param1=usestats{i};
        Param2=usestats{j};
        mystat=cellfun(@(a) mean(a.(Param1),'omitnan'),useSession.CallStats(:));
        mystat2=cellfun(@(a) mean(a.(Param2),'omitnan'),useSession.CallStats(:));
        figure;
        sp=subplot(1,2,1);
        scatter(mystat,mystat2,5,daycolor,'filled');
        xlabel(Param1); ylabel(Param2);
        cb=colorbar('Ticks',linspace(0,1,8),'TickLabels',4:2:20); cb.Label.String='Day';
        sp(2)=subplot(1,2,2);
        scatter(mystat,mystat2,5,genocolor,'filled');
        xlabel(Param1); ylabel(Param2);
        cb=colorbar('TickLabels',genonames,'Ticks',0.002:0.004:0.014);
        cb.Label.String='Day'; cb.Limits=[0,0.0156];
        colormap(gca,'lines')
        linkaxes(sp);
    end
end
 % do this but just plot males to see
        
        
        
mystat=cellfun(@(a) mean(a.(Param1),'omitnan'),useSession.CallStats(:));
mystat2=cellfun(@(a) mean(a.(Param2),'omitnan'),useSession.CallStats(:));
figure;
sp=subplot(1,2,1);
scatter(mystat,mystat2,5,daycolor,'filled');
xlabel(Param1); ylabel(Param2);
cb=colorbar('Ticks',linspace(0,1,8),'TickLabels',4:2:20); cb.Label.String='Day';
sp(2)=subplot(1,2,2);
scatter(mystat,mystat2,5,genocolor,'filled');
xlabel(Param1); ylabel(Param2);
cb=colorbar('TickLabels',genonames,'Ticks',0.002:0.004:0.014);
cb.Label.String='Day'; cb.Limits=[0,0.0156];
colormap(gca,'lines')
linkaxes(sp);



%%


%% input data for the clustering algorithms
%{
First, lets understand the unput data that DeepSqueak uses
1. Entropy: geomean/mean of the full window
2. ridgeTime: the indices in the call box of the ridge (calced by either
	bs power or entropy (I use blobs instead)
3. ridgeFreq: they use a single frequency, i may be able to get away with a
    few
4. smoothed ridge frequency (using a rlowess smoothing algo, using .1 input
5. a filtered image using imgradientxy on the raw image I
6. signal to noise: mean 1-entropy of the ridge
7. bein time (of the ridge)
8. end time (real time)
9. DeltaTime= duration of ridgeline
10. Low Freq- of ridge
11. High Freq of ridge
12. DeltaFreq of above two
13. stdev= std freqscale*smoothed ridgefreq (scale must be the df)
14. slope is some odd thing... line 116... i think its a matrix version of
    the contour slope
15. MaxPower= mean ridgepower
16. Power=ridgepower (the vector)
17. RidgeFreq=the actual freqs of the ridge
18. PeakFreq= freq at highest power of ridge (good idea!)
19. Sinuosity= line distance of the ridgeline over the duration of the
    ridgeline ** good one too, i bet we can weighted average across all
    ridgelines

FROM THIS: they only use a few stats for each clustering algorithm

The output:
1. Spectrogram
2. lower freq of box
3. delta time
4. xfreq (Freq pts)
5. xTime (time points)
6. file path
7. perFileCallID
8. power
9. box(4) which is height

For Kmeans:
1. reshapedX
2. zscore slope of 

%}

%% this is the beginning of the Procedure file
%
% produced by zach and jay
%
%% for zachs data heres how they got analyzed:

% zach has a file that he loads first... but he also exported the excel
% data

% this imports the excel calls into a matlab variable
edit importdeepsqueakexcel
% this parses the filenames
edit preprocess_deepsqueak
% this is inside the above to pull data
edit scrape_fileinfo
%%
%%%%%%%%  now analyze data %%%%%%%%%%%%%
edit calls_count_analysis
edit full_plot_analysis