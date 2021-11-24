%{
This is the analysis sketchpad of the usv data from start to finish
 I will start with blocks of sketch code here and then offload them into
 functions or other scripts if they become useful.

the preprocessing pipeline was as follows:

1. DeepSqueak:
-data were acquired, see acquisition protocol
-a systematically random selection of files were pulled and manually
curated for training data.
-a network was built on those calls using DeepSqueak v 3.0.1
-the full detection algorithm was then performed on the full gambit of data
using both the rat detector yolo r1 and rat detector JHB networks together
-a home made denoiser was then applied to ALL FILES

2. using USVseg, the problem basically is that they fucked up the
multitaper spectrogram, so it doesnt work... its possible that we could use
chronux bit i havent tried.  I am rather convinced that multitaper isnt
really useful unless you have a lot of noise, so im willing ot kill it for
now
-so I would hack off USVSeg into my own pipeline, run all the files,
instead of spitting out a csv, i can add them in their own table to the
struct, and go from there.

The idea would be to okay all usvs that get caught by both, maybe using the
fattest window first, and then those that only get caught by one, i'd have
to peek and see what tehy miss.

3.vocalMat:
   



Manual edits:
First some observations:

-the bounding boxes tend to miss the start and end of the calls, leaving
loud bits outside the box

- The contour isnt quite right sometimes, need to check how they threshold
to get a contour.
    The contour math uf deepsqueak:
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
* the scattoni paper has some ok math on how to identify different usvs.
I'll have to do some work on this though.

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

% this basically takes the deepsqueak data, compares them to usvseg data,
% and then removes calls that suck...
edit usvsegDetect

% and to implement the usvsegDetect:

wavfiledir='G:\USV data\Raw\wav';
wavfilenames=getAllFiles(wavfiledir);
fileinfo=dir(wavfiledir); fileinfo([fileinfo.isdir])=[];
outdir='G:\USV data\segData';
if ~isfolder(outdir)
    mkdir('G:\USV data\segData');
end
nameadd='_callData';
bigclock=tic;
wb=waitbar(0,'Starting to build each session');
for i=427:length(wavfilenames)
    smallclock=tic;
    audiodata=audioinfo(wavfilenames{i});
    
    sessdata=rmfield(fileinfo(i),{'bytes','isdir'});
    myfilename=sessdata.name;
    unders=find(myfilename=='_' | myfilename=='-');
    sessdata.Cohort=str2double(myfilename(2:unders(1)-1)); % cohort
    sessdata.Age=str2double(myfilename(unders(1)+2:unders(2)-1)); % age
    % pull rat name from file name
    sessdata.Ratname=myfilename(unders(2)+1:end-4);
    
    [spect,thrshd,params,onoffset,onoffsetm,blobs]=usvsegDetect(audiodata);
    % spect is the flattened centered im, thrshd is the pix that pass
    % threshold after masking
    
    segCalls=table(onoffset(:,1),onoffset(:,2),ones(size(onoffset,1),1),...
        'VariableNames',{'onsetTime','offsetTime','Accept'});
    
    calldata=getCallStats(segCalls,params,spect,thrshd);
    
    % reduce the data size of the spect
    spect=uint8(rescale(spect)*256);
    save(fullfile(outdir,[sessdata.name(1:end-4) nameadd]),...
        'audiodata','params','spect','blobs','segCalls','calldata');
    waitbar(i/length(wavfilenames),wb,sprintf('Its taken %d minutes, Likely %d more',...
        round(toc(bigclock)/60),round(toc(smallclock)*(length(wavfilenames)-i)/60)));
end

% the next step i think would be to get all the binary images and turn them
% into a smaller vector (like they do in deepsqueak using a vae)

%% now to drop these into smaller dimensions using a VAE

% i will generate a new latentscore variable for each day that is the vae
% deconstruction of the data.  I'll have to save the vae encoder and
% decoder data so that I can use it to reconstruct centroids and examples.

% first generate the input dataset.
% this time try the black and white images
inputfolder=uigetdir('Select folder of usvSeg Data');
statfiles=dir(inputfolder);
statfiles=statfiles(~[statfiles.isdir]);
load('USVmetadataC1-9.mat');
USVrecordings=mergevars(USVrecordings,{'Animal','EarMarkings','LegMarkings','ToeMarkings'},'NewVariableName','Aliases');

USVlegend=table;
for i=1:length(statfiles)
    filedata=statfiles(i);
    % attach metadata to each stat file;
    sessdata=rmfield(filedata,{'bytes','isdir'});
    myfilename=sessdata.name;
    unders=find(myfilename=='_' | myfilename=='-');
    sessdata.cohort=str2double(myfilename(2:unders(1)-1)); % cohort
    sessdata.age=str2double(myfilename(unders(1)+2:unders(2)-1)); % age
    % pull rat name from file name
    sessdata.ratname=string(myfilename(unders(2)+1:unders(3)-1));
    sessdata.name=string(sessdata.name);
    % now match my animal
    sameCohort=USVrecordings.Cohort==sessdata.cohort;
    sameName=contains([USVrecordings.Aliases],sessdata.ratname,'IgnoreCase',true);
    % find this rat, and then add its info
    sameRat=find(sum([sameCohort sameName],2)>1);
    if length(sameRat)==1
        filedata=[struct2table(sessdata) USVrecordings(sameRat,:)];
    elseif length(sameRat)>1
       fprintf('Ran into multiple matches for %d %s \n', sessdata.cohort,sessdata.ratname);
       fprintf('   %s \n',USVrecordings.Aliases(sameRat,:)');
       fprintf('Choosing first \n \n')
       filedata=[struct2table(sessdata) USVrecordings(sameRat(1),:)];
    else
       fprintf('Ran into no matches for %d %s \n', sessdata.cohort,sessdata.ratname);
    end
    USVlegend(i,:)=filedata;
end

%%

% lets take just the males, and just between P4 to P10
% this is hardcoded for now, but suffice to say its going to be 600 msec
% images from 15 khz to 90 khz
runSess=USVlegend(lower(USVlegend.Sex)=='m' & USVlegend.age<13,:);
allcalls=[]; allCallinfo=[];
wb=waitbar(0,'concatenating ALL the images');
% what percentage of images will we need to classify these guys????
% probably like half
bigclock=tic;
for i=1:height(runSess)
    smallclock=tic;
    load(fullfile(runSess.folder(i,:),runSess.name(i)),'blobs','segCalls','params');
    okfreqs=params.fvec>15000 & params.fvec<90000;
    badcalls=segCalls.onsetTime<5 | segCalls.offsetTime>min([180 size(blobs,2)*params.timestep])-5;
    segCalls=segCalls(~badcalls,:);
    mycalls=false(128,128,1,height(segCalls));
    callcenters=mean(table2array(segCalls(:,1:2)),2);
    centerinds=interp1(params.tvec,1:length(params.tvec),callcenters,'nearest');
    calldurs=(segCalls.offsetTime-segCalls.onsetTime)/params.timestep;
    
    % go forward 200 msec, and forward 200 msec, then resize to 256
    for j=1:length(callcenters)
    
        % were going from many to 128x128 x 1 images
        mycalls(:,:,:,j)=imresize(blobs(okfreqs,(centerinds(j)-round(calldurs(j)/1.9)):(centerinds(j)+round(calldurs(j)/1.9))),[128,128]);
    end
    allCallinfo=[allCallinfo; [calldurs ones(size(calldurs,1),1)*i (1:length(calldurs))']]; % add call duration, then the sessnum and callnum
    allcalls=cat(4,allcalls,mycalls);
    waitbar(i/height(runSess),wb,sprintf('Has taken %d mins, will prolly go %d more mins',...
    round(toc(bigclock)/60),round(toc(smallclock)/60*(height(runSess)-i))));
end
close(wb);

%%
% use create_tsne callback to generate this code

edit generate_VAE_encoder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% USE [encoderNet,decoderNet,data] = generate_VAE_encoder(imageStack)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now run the encoder decoder...
% Load the network model

% then yu can run tsne or umap on this and cluster
% first will need to pull a random set of usvs to generate our model.  This
% is throttled by the size of the gpu i guess...

batchsize=20000; 

% 1. generate random sampling
randpull=randperm(size(allcalls,4));
trainingstack=allcalls(:,:,:,randpull(1:batchsize));

% 2. build encoder and decoder using those images
[encoderNet,decoderNet,data] = generate_VAE_encoder(trainingstack);

clear trainingstack;
% 3. run the encoder throguh ALL images
batchinds=[[1 batchsize+1:batchsize:size(allcalls,4)]' [batchsize:batchsize:size(allcalls,4) size(allcalls,4)]'];
callData=[];
for i=1:size(batchinds,1)
    [~,~,batchData]=generate_VAE_encoder(allcalls(:,:,:,batchinds(i,1):batchinds(i,2)),encoderNet,decoderNet,1);
    callData=[callData; batchData];
end

% 4. put the callduration back in
% 32 dims, 33rd is callduration, and allCallinfo is session and call number
% from runSess;
callData=[callData allCallinfo(:,1)];

% 5. see if any of these predict genotype
% gather our genotypes
for i=1:max(allCallinfo(:,2))
    allCallinfo(allCallinfo(:,2)==i,4)=contains(runSess.Genotype(i),'fx');
end    


[b]=fitglm(callData,allCallinfo(:,4));

% 6. if not, cluseter calls, see if one genotype dominates one cluster, or
% whether one call varies by genotype
mycolors=lines(2);
for i=2: size(callData,2)
    for j=1:i-1
        figure;
        scatter(callData(:,i),callData(:,j+6),4,mycolors(allCallinfo(:,4)+1,:),'filled');
    end
end



%% try tsne on these animals
% some very basic questions:
%{
1. do these calls cluster at all in any way
2. do these calls stratify by.... animal, genotype, time....


%}

%% or umap
if ~contains(path,'C:\Users\Jadhavlab\Documents\gitRepos\USV_Analysis\Umapfx\umap')
    addpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\USV_Analysis\Umapfx\umap'));
end
[embed2,umap2,clusterid] = run_umap(data);
rmpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\USV_Analysis\Umapfx\umap'));
%%
ratNames=string();
for i=1:length(allCallinfo)
    ratNames(i)=runSess.Aliases(allCallinfo(i,2),1);
end
[a,b,c]=unique(ratNames);
mycolors=parula(length(a));
for perplexity=1:8
embed = tsne(callData,'Verbose',1,'Perplexity',perplexity*4+16);

embed = (embed - min(embed)) ./ (max(embed)-min(embed));
embedY = 1-embed(:,2); % flip Y coordinates so the images looks like the UMAP figure
embedX = embed(:,1);

figure;
scatter(embedY,embedX,5,mycolors(c,:),'filled');
title(sprintf('perplexity %d', perplexity*6));
end

%% using deepsqueak...
callfiledir='G:\USV data\Detections';
if ~isfolder(callfiledir)
    callfiledir=uigetdir;
end
callfiles=getAllFiles(callfiledir);
wb=waitbar(0,'starting');
tic;
% this will take roughly 12 hours to run, or overnight
for i=1:length(callfiles)
    thisf=tic;
    load(callfiles{i});
    
    % this has variables: CALLS AUDIODATA DETECTION_METADATA SEGCALLS
    % we probably want to now save the spectrogram as well as the bwimage
    % (but only for valid blobs)
    
    % may be able to optimize with guparray
    [spect,thrshd,params,onoffset,onoffsetm,blobs]=usvsegDetect(audiodata);
    % spect is the flattened centered im, thrshd is the pix that pass
    % threshold after masking

    segCalls=table(onoffset(:,1),onoffset(:,2),ones(length(onoffset),1),...
        'VariableNames',{'onsetTime','offsetTime','Accept'});
    CallStats=getCallStats(segCalls,params,spect,thrshd);

    %{
    stream = RandStream('mlfg6331_64');  % Random number stream
    options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
    % need to do this with different number of means and gather some stats
    % on cluster quality (prob iso and l-ratio)
    [idx,C,sumd,D] = kmeans(table2array(allCalls),5,'Options',options,'MaxIter',10000,...
        'Display','final','Replicates',10);
    %}
    

    
    %{
    for i=1:max(idx)
        figure;
        [inds]=find(idx==i);
        toplot=randperm(length(inds));
        for j=1:min([length(toplot) 9])
            subplot(3,3,j)
            mybounds=[Calls.realWin(inds(toplot(j)),1) Calls.realWin(inds(toplot(j)),2)];
            onsetI=max([1 find(params.tvec<=mybounds(1),1,'last')]);
            offsetI=min([find(params.tvec>mybounds(2),1,'first') length(params.tvec)]);
            imagesc((onsetI:offsetI)*tstep,params.fvec,thrshd(:,onsetI:offsetI))
        end
        sgtitle(sprintf('Category %d',i));
    end
    allCalls2=allCalls; allCalls2.specGram=[];
    [coeff,score,latent] = pca(zscore(table2array(allCalls2)));
    figure; scatter3(score(:,1),score(:,2),score(:,3));
    %}

    save(callfiles{i},'Calls','audiodata','params','spect','blobs','segCalls','CallStats');
    fprintf('took %.f seconds\n',toc(thisf));
    waitbar(i/length(callfiles),wb,...
        sprintf('runing now, probably %.2f mins left',toc/i/60*(length(callfiles)-i)));
end
close(wb);


%%

edit removeSoftCalls

try
    load('USVmetadata2021-10-15.mat');
catch
    edit USVmetaAggregator 
end

edit mergeMetaCalls

edit ClusterCalls

%% just curious about weight
weightmat={}; allweights={};
for i=1:5 %max(USVSession.Cohort)
    uniquedays=unique(USVSession.Age(USVSession.Cohort==i));
    allweights{i}=uniquedays;
    for j=1:length(uniquedays)
        weightmat{i}(:,j)=USVSession.Weight(USVSession.Cohort==i & USVSession.Age==uniquedays(j));
    end
end
figure;
temp={1,2,3,4,5,6}; cohortnames=cellfun(@(a) {['Cohort ' num2str(a)]}, temp);
clear eb;
for i=1:length(weightmat)
    eb(i)=errorbar(allweights{i},nanmean(weightmat{i}),nanstd(weightmat{i}));
    hold on;
end
xlabel('Age pnd'); ylabel('Weight, gm')
legend(eb,cohortnames)
% it appears that cohort 5 is older than reported...
%%

load('C:\Users\John Bladon\Desktop\USVdataFULL2021-10-15.mat');
% this is what deepsqueak gives me to work with...

%% now that we have a datatable, lets see what falls out
% first lets plot everything by sex/genotype and time
oktouse=cellfun(@(a) ~isempty(a), USVSession.CallStats(:)) & ...
    USVSession.Cohort~=6;
useSession=USVSession(oktouse,:);
allstats=useSession.CallStats{1}.Properties.VariableNames;
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

figure;
clear eb
for gt=1:4 % 4 genotypes
    daymean=accumarray(useSession.Age(idx==gt),useSession.Weight(idx==gt),[],@nanmean,nan);
    daydev=accumarray(useSession.Age(idx==gt),useSession.Weight(idx==gt),[],@nanstd,nan);
    alldays=unique(useSession.Age(idx==gt));
    eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
        'Color',gtcolors(gt,:));
    hold on;
end
 xlim([3 21]);
legend(eb,genonames);
ylabel('Weight, gms'); xlabel('Age, pnd');

% lets see if this is real
%[p,tbl,stats]=anovan(useSession.Weight,{useSession.Age, idx},'continuous',1,...
%    'varnames',{'Age','Genotype'},'model','full');

[p,tbl,stats]=anovan(useSession.Weight,{useSession.Age, idx==3},'continuous',1,...
    'varnames',{'Age','Genotype'});


%% so I think we should plot the p(call) for each call type for each animal

allcalls=cellfun(@(a) a.Label, useSession.CallStats(:), 'UniformOutput', false);
everycall=cell2mat(cellfun(@(a) str2double(a), allcalls, 'UniformOutput', false));
everycall(isnan(everycall))=26; % there are 26 types of call (26 is unidentified)
uniquecalls=unique(everycall);
savedir='E:\Brandeis datasets\FMR1 Project Data\USV data\Figs 2021-10-18';
goodcalls=[2 4 5 6 7 8 9 16 17 18];
for i=1:length(goodcalls)
    figure;
    mystat=cellfun(@(a) mean(str2double(a.Label)==goodcalls(i),'omitnan'),useSession.CallStats(:));
    %sp=subplot(3,4,i+1);
    for gt=1:4 % 4 genotypes
        daymean=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@nanmean,nan);
        daydev=accumarray(useSession.Age(idx==gt),mystat(idx==gt),[],@SEM,nan);
        alldays=find(~isnan(daymean));
        eb(gt)=errorbar(alldays+rand*.5-.25,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
            'Color',gtcolors(gt,:));
        hold on;
    end
    % pull out random images for those clusters
    %sp(2)=
    xlabel('Postnatal Age'); ylabel('likelihood of call ID#');
    box off; title(sprintf('call %d',goodcalls(i)));
    legend(eb,genonames); xlim([3 21]);
    %saveas(gcf,fullfile(savedir,sprintf('callnumber_%d',goodcalls(i))),'png');
    %close(gcf);
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
burstprob=[];
wb=waitbar(0,'starting');
for i=1:height(useSession)
    allISI=diff([useSession.CallStats{i}.EndTimes(1:end-1) useSession.CallStats{i}.BeginTimes(2:end)],1,2);
    [y,x]=histcounts(allISI,0:.005:.4);
    [~,peakind]=max(y); isiPeak(i)=xcenters(peakind);
    [y,x]=histcounts(useSession.CallStats{i}.CallLengths,xedges);
    [~,peakind]=max(y); callDurPeak(i)=xcenters(peakind);
    burstprob(i)=nanmean(allISI<.2);
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
figure;
%sp=subplot(3,4,i+1);
for gt=1:4 % 4 genotypes
    daymean=accumarray(useSession.Age(idx==gt),burstprob(idx==gt),[],@nanmean,nan);
    daydev=accumarray(useSession.Age(idx==gt),burstprob(idx==gt),[],@nanstd,nan);
    alldays=find(~isnan(daymean));
    eb(gt)=errorbar(alldays,daymean(~isnan(daymean)),daydev(~isnan(daydev)),...
        'Color',gtcolors(gt,:));
    hold on;
end
xlabel('Postnatal Age');
ylabel('Call Bout Probability, seconds')
legend(eb,genonames)
% 

%% cluster early day animals, cluster late day animals





%% plot certain characteristics against eachother
useSession=USVSession(lower(USVSession.Sex)=='m'& ...
    cellfun(@(a) ~isempty(a), USVSession.CallStats(:)) ,:);
[alldays,~,dayind]=unique(useSession.Age);
daycolor=parula(length(alldays));
daycolor=daycolor(dayind,:);

[genotypes,~,idx]=unique(useSession.Genotype);
genonames={'fx M','wt M'};
gtcolors=lines(2);
genocolor=gtcolors(idx,:);




usestats=allstats(7:17); %just use all these ones...
for i=1:length(usestats)-1
    for j=i+1:length(usestats)
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
        cb=colorbar('TickLabels',genonames(1:2),'Ticks',[.002 .006]);
        cb.Label.String='Genotype'; cb.Limits=[0,0.007];
        colormap(gca,'lines')
        linkaxes(sp);
    end
end
 % do this but just plot males to see

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