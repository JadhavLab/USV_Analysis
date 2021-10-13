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

%% reject calls that are too soft

% first import a tool file
load('G:\USV data\Detections\C1_P6_1BL 2021-10-07  1_56 PM.mat');
% variables are audiodata, Calls, detection_metadata

% audiodata is how to get the data...
% so we need to whiten the whole dataset.  to do this, we zscore each freq
% independently.  Luckily i have a gpu that can handle that...

audiodata.Data=gpuArray(audioread(audiodata.Filename)); % pull all samples

% now load the full spectrogram onto gpu
audiodata.windowsize = round(audiodata.SampleRate*0.0032); % hardcoded in for large spectrogram at least
audiodata.noverlap = round(audiodata.SampleRate*.0016);
audiodata.nfft = round(audiodata.SampleRate*.0032);
%% Make the spectrogram
span=1:10000;



[spect.raw, f, t] = spectrogram(audiodata.Data,audiodata.windowsize,audiodata.noverlap,audiodata.nfft,audiodata.SampleRate,'yaxis');
% get the abs value for amplitude
spect.samp=log(abs(spect.raw));
% adjust an image
spectmedian=mean(spect.samp,2);
spectspread=std(spect.samp,[],2);
spect.squeakim=SmoothMat2(spect.samp(:,span)-spectmedian,[10 10],2);

% now get median and STD
figure; 
sp=subplot(2,2,1);
imagesc(t(span),f,spect.squeakim);
set(gca,'Ydir','normal','CLim',[-1 2],'YLim',[20000,100000]);
% establish a threshold
subplot(2,2,2);
plot(f,spectmedian);
yyaxis right; semilogy(f,spectspread);
set(gca,'Xlim',[20000,100000])

% now threshold
sp(2)=subplot(2,2,3);
imagesc(t(span),f,(spect.samp(:,span)>spectmedian+spectspread*2))
set(gca,'Ydir','normal','YLim',[20000,100000]);


subplot(2,2,4);
ecdf(linearize(abs(spect.raw(1:find(f<100000,1,'last'),1:10000))))
set(gca,'XScale','log');
hold on; plot(repmat(mean(spectmedian+2*mean(spectspread)),1,2),[0 1]);

linkaxes(sp,'x','y');

%%

% to put this into practice:

% 1. for each file, get the full spect
% 2. calc the median and spread
% 3. go through each call and define the blobs that are large and exceed
% the threshold
% 4. if there are none, or they are too small, remove the call.

%filename='G:\USV data\Detections\C1_P6_1BL 2021-10-07  1_56 PM.mat';

mydir='G:\USV data\Detections'; %uigetdir([],'Get folder with all the raw call mats');

allfiles=getAllFiles(mydir,'.mat');
%%
wb=waitbar(0,'Starting up the curation process');
for fi=420:length(allfiles)

    load(allfiles{fi});
    % only go to 3 minutes
    % first only grab calls that start within 180 seconds
    Calls=Calls(Calls.Box(:,1)<180,:);
    
    if ~ismember('blobs',Calls.Properties.VariableNames)
        
        gpuDevice(1); % reset the gpu by clearing it
        audiodata.Data=gpuArray(audioread(audiodata.Filename,[1 min([180 audiodata.Duration])*audiodata.SampleRate])); % pull all samples
        
        % now load the full spectrogram onto gpu
        audiodata.windowsize = round(audiodata.SampleRate*0.0032); % hardcoded in for large spectrogram at least
        audiodata.noverlap = round(audiodata.SampleRate*.0016);
        audiodata.nfft = round(audiodata.SampleRate*.0032);
        spect=struct('raw',[]);
        [spect.raw, f, t] = spectrogram(audiodata.Data,audiodata.windowsize,audiodata.noverlap,audiodata.nfft,audiodata.SampleRate,'yaxis');
        
        % clip what we dont want to see, anything below 10 khz or above 100 khz
        freqbounds=double(find(f>10000,1,'first'):find(f<100000,1,'last'));
        spect.samp=log(abs(spect.raw(freqbounds,:)));
        spect.f=f(freqbounds);
        spectmedian=gather(mean(spect.samp,2));
        spectspread=gather(std(spect.samp,[],2));
        spect=rmfield(spect,'raw');
        t=gather(t);
        
        % load the im for each call
        for i=1:height(Calls)
            % start 50 ms before call (USVseg has a max of a 30 msec pause for it
            % to be the same call, any gap larger and its a different call
            boxstart=max([1 find((Calls.Box(i,1)-.05)>t,1,'last')]);
            boxend=min([length(t) find((Calls.Box(i,3)+.05+Calls.Box(i,1))<t,1,'first')]);
            
            % if call goes past start or end of recording just reject
            if boxstart==1 || boxend==length(t)
                Calls.Accept(i)=0;
                continue
            end
            
            callspect=spect.samp(:,boxstart:boxend);
            callT=t(boxstart:boxend);
            % now gather blobs
            normalizedspect=(callspect-spectmedian)./spectspread;
            
            
            
            callblobs=regionprops(gather(normalizedspect>2),gather(normalizedspect),...
                'PixelIdxList','Area','BoundingBox','MeanIntensity',...
                'Centroid','Orientation','Eccentricity','Extent','EulerNumber','FilledImage');
            % has to have at least
            % blobs have to be:
            % at least 6 msec (4 pix horizontal),
            %at least 15 pixels total,
            % mean intensity must be above threshold of 2.2 sd above mean
            % orientation is greater than say 80*, or basically vertical
            % must not be 'holey' e.g. filled image cant be more than 110% of
            % holes image
            % most have a centroid higher than 180kHz
            % blob must not be through start or end of the call box
            
            okidx= cellfun(@(a) a(3)>4, {callblobs.BoundingBox}) &...
                [callblobs.Area]>15 & ...
                [callblobs.MeanIntensity]>2.2 & ...
                ~([callblobs.Orientation]>80 | [callblobs.Orientation]<-80)&...
                (cellfun(@(a) sum(a(:)), {callblobs.FilledImage})./[callblobs.Area]<1.1) & ...
                cellfun(@(a) spect.f(round(a(2)))>18000, {callblobs.Centroid}) &...
                cellfun(@(a) a(1)>5 && (a(1)+a(3))<(size(callspect,2)-5),{callblobs.BoundingBox});
            callblobs=callblobs(okidx);
            
            verbose=0;
            if verbose
                mf=figure('Position',[1500,420,1000,920]); subplot(2,2,[1 3]);
                imagesc(callT,spect.f,normalizedspect);
                subplot(2,2,2);
                imagesc(callT,spect.f,(normalizedspect>2)+(normalizedspect>2.2));
                blobIM=zeros(size(normalizedspect,1),size(normalizedspect,2));
                for bl=1:length(callblobs)
                    blobIM(callblobs(bl).PixelIdxList)=bl;
                end
                subplot(2,2,4);
                imagesc(callT,spect.f,blobIM);
                keyboard
                close(mf);
            end
            if isempty(callblobs)
                Calls.Accept(i)=0;
                continue
            else
                Calls.blobs{i}=callblobs;
            end
        end
    end
% now save out
save(allfiles{fi}, 'Calls', '-append');
waitbar(fi/length(allfiles),wb,sprintf('Running Files Now! file %d of %d',fi,length(allfiles)));
end
waitbar(1,wb,'All Done! Kill this box');
       
%%
    
    if verbose
        figure; subplot(1,2,1);
        imagesc(callT,spect.f,normalizedspect);
        subplot(1,2,2);
        imagesc(callT,spect.f,(normalizedspect>1.8)+(normalizedspect>2.2));
    end
%These are all usvseg techniques for getting decent thresholded data.  I
%find these are not that useful and do not demostrably change the
%spectrogram
%{

% liftering the data

figure;
sp=subplot(2,1,1);
imagesc(log(abs(spect.raw(:,span))));
set(gca,'ydir','normal','Clim',[-5 1]);
liftercutoff = 3; % fixed parameter
fftsize = (size(spect.raw(:,span),1)-1)*2;
cep = fft([spect.raw(:,span);flipud(spect.raw(2:end-1,span))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);
sp(2)=subplot(2,1,2);
imagesc(log(abs(liftered(:,span))));
set(gca,'ydir','normal','Clim',[-5 1]);
linkaxes(sp);
% liftering seems to do nothing demonstrable to the data

% Median centering


% median centering
med = median(liftered,2);
liftmed = liftered-repmat(med,1,size(liftered,2));
figure;
sp=subplot(2,1,1);
imagesc(log(abs(spect.raw(:,span))));
set(gca,'ydir','normal','Clim',[-5 1]);
title('raw');
sp(2)=subplot(2,1,2);
imagesc(log(abs(liftmed(:,span))));
set(gca,'ydir','normal','Clim',[-5 1]);
linkaxes(sp);
title('liftered median centered');
% potentially improved the data may have lowered the noise

% movmedian smoothed

if exist('movmedian')==2
    fltnd = movmedian(liftmed,5);
else % for R2015b or earlier
    pad1 = repmat(liftmed(1,:),2,1);
    pad2 = repmat(liftmed(end,:),2,1);
    padded = [pad1;liftmed;pad2];
    fltnd = zeros(size(liftmed));
    for n=1:size(liftmed,1)
        fltnd(n,:) = median(padded(n+(0:4),:));
    end
end

figure;
sp=subplot(2,1,1);
imagesc(log(abs(spect.raw(:,span))));
set(gca,'ydir','normal','Clim',[-5 1]);
title('raw');
sp(2)=subplot(2,1,2);
imagesc(log(abs(spect.raw(:,span)-fltnd*2)));
set(gca,'ydir','normal','Clim',[-5 1]);
linkaxes(sp);
title('liftered median centered median smoothed');
%}

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
8.end time
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