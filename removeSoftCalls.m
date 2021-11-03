function [Calls,segCalls,squeakCalls]=removeSoftCalls(Calls,onoffsetm,onoffset,params,thrshd,spect)
% [Calls,segCalls,squeakCalls]=removeSoftCalls(Calls,onoffsetm,onoffset,params,thrshd,spect)
% function removeSoftCalls compares calls from DeepSqueak and usvSeg and
% determines which are not real calls.
% INPUTS:
%   Calls
%   onoffsetm
%   onoffset
%   params
%   thrshd
%   spect
% OUTPUTS:
%   Calls
%   segCalls
%   squeakCalls

% Algorithm:
%   1. merge all calls by window; use the fattest window for all calls
%   2. review all calls that are from a single detector algorithm
%   3. reconstitute the spectrogram there, and assay if there is enough
%       blobs of intensity to allow the call
%   4. Output new combined Call timestamps.

% JHB 11-3-2021

squeakCalls=Calls;
segCalls=[onoffsetm onoffset ones(length(onoffsetm))];
squeakTimes=[squeakCalls.Box(:,1) squeakCalls.Box(:,1)+squeakCalls.Box(:,3)]; squeakTimes(:,3)=1;
segTimes=onoffsetm; segTimes(:,3)=2;
allCalls=sortrows([squeakTimes; segTimes],1);
allCalls(:,2)=sort(allCalls(:,2)); % independently sort ends
% syllables whose margins are overlapped are integrated in onoffsetm but not in onoffset
overlap = find((allCalls(2:end,1)-allCalls(1:end-1,2))>=0); % use only those onto which the offset is before or at next onset
onsetTime = [allCalls(1);allCalls(overlap+1,1)]; % if they are, notch the on and offs out of the center
offsetTime = [allCalls(overlap,2); allCalls(end,2)];
callID= [[allCalls(1,3); allCalls(overlap+1,3)] [allCalls(overlap,3); allCalls(end,3)]];

Calls=table(onsetTime,offsetTime,callID,ones(length(onsetTime),1),...
    'VariableNames',{'onsetTime','offsetTime','squeak1seg2','Accept'});
% so this is all the calls without copies, now we go through each and
% validate? how would we validate though... first lets see how the
% algos do.. looks like 13% in this file are exclusive to one algo
% 6% are deepsqueak, and 7% are usvseg.  I think we can consider 87%
% validated, and then go after the 6 and 7%

orphanID=find(diff(callID,1,2)==0);
Calls.realWin(1,:)=[nan nan]; % initiate realwin
for id=1:length(orphanID)

    % sanity checking here


    onsetI=max([1 find(params.tvec<=onsetTime(orphanID(id))-.03,1,'last')]);
    offsetI=min([find(params.tvec>offsetTime(orphanID(id))+.03,1,'first') length(params.tvec)]);

   
    % okay this is a good example of noise, so lets use blob analysis
    maskbounds=[onsetTime(orphanID(id)), offsetTime(orphanID(id)),...
        15000, 90000];
    mask=ones(size(spect,1),offsetI-onsetI);
    mask(:,params.tvec(onsetI:offsetI)<maskbounds(1) | params.tvec(onsetI:offsetI)>maskbounds(2))=0;
    mask(params.fvec<maskbounds(3) | params.fvec>maskbounds(4),:)=0;


    callblobs=regionprops(logical(thrshd(:,onsetI:offsetI)),spect(:,onsetI:offsetI),...
        'PixelIdxList','Area','BoundingBox','MeanIntensity',...
        'Centroid','Orientation','Eccentricity','Extent','EulerNumber','FilledImage');
    % has to have at least
    % blobs have to be:
    % at least 2 msec, should be about 5 pix across
    durmin=round(.002/params.timestep);
    % mean intensity must be above threshold of 2.2 sd above mean
    % orientation is greater than say 80*, or basically vertical
    % must not be 'holey' e.g. filled image cant be more than 110% of
    % holes image
    % most have a centroid higher than 180kHz
    % blob must not be through start or end of the call box

    okidx= cellfun(@(a) a(3)>durmin, {callblobs.BoundingBox}) &...
        [callblobs.MeanIntensity]>2.2 & ...
        ~([callblobs.Orientation]>80 | [callblobs.Orientation]<-80)&...
        (cellfun(@(a) sum(a(:)), {callblobs.FilledImage})./[callblobs.Area]<1.1) & ...
        cellfun(@(a) a(1)>5 && (a(1)+a(3))<(offsetI-onsetI-5),{callblobs.BoundingBox});
    callblobs=callblobs(okidx);
    
    % reconstitute the blobs only:
    callImage=zeros(size(thrshd(:,onsetI:offsetI)));
    for bl=1:length(callblobs)
        inBox=ceil(callblobs(bl).BoundingBox);
        callImage(inBox(2):inBox(2)+inBox(4)-1,inBox(1):inBox(1)+inBox(3)-1)=...
            callblobs(bl).FilledImage;
    end
    
     % image the thing;
    figure; sp=subplot(4,1,1);
    imagesc(params.tvec(onsetI:offsetI),params.fvec,spect(:,onsetI:offsetI));
    sp(2)=subplot(4,1,2);
    imagesc(params.tvec(onsetI:offsetI),params.fvec,thrshd(:,onsetI:offsetI));
    sp(3)=subplot(4,1,3);
    imagesc(params.tvec(onsetI:offsetI),params.fvec,thrshd(:,onsetI:offsetI).*mask);
    sp(4)=subplot(4,1,4);
    imagesc(params.tvec(onsetI:offsetI),params.fvec,callImage);



    %if isempty(callblobs)
    if isempty(callblobs)
        Calls.Accept(orphanID(id))=0;
        Calls.realWin(orphanID(id),:)=[nan nan];
    else
        keyboard;
        Calls.Accept(orphanID(id))=1;
        % get the start and end of all the okay blobs
        blobinds=cell2mat({callblobs.BoundingBox}');
        Calls.realWin(orphanID(id),:)=[params.tvec(onsetI+min(ceil(blobinds(:,1)))) ...
            params.tvec(onsetI+max(ceil(blobinds(:,1))+blobinds(:,3)-1))];
    end

end
end


%% removeSoftCalls old code

%%%% reject calls that are too soft %%%%%%%%%
%{
% first import a tool file
load('G:\USV data\Detections\C3-P20-T8 2021-10-07  3_26 PM.mat');
% variables are audiodata, Calls, detection_metadata

% audiodata is how to get the data...
% so we need to whiten the whole dataset.  to do this, we zscore each freq
% independently.  Luckily i have a gpu that can handle that...

audiodata.Data=gpuArray(audioread(audiodata.Filename)); % pull all samples

% compare the spectrogram settings here:


% now load the full spectrogram onto gpu
% this can go much smaller, i use .8 msec bins, could prolly do .6
audiodata.windowsize = round(audiodata.SampleRate*0.0032); % hardcoded in for large spectrogram at least
audiodata.noverlap = round(audiodata.SampleRate*.0016);
audiodata.nfft = round(audiodata.SampleRate*.0032); % 1229 but it will go up to too high hz




winsize=.0006;
winbins=round(audiodata.SampleRate*winsize);

stepsize=round(winbins/2); % use a 2x overlap

nfft=1229*2; % use double, you kill the top half anyways
% or use fvec....
fvec=0:120:100000;
%}

%% this section is the old thresholding that the current method is based on
% see below:

% to put this into practice:

% 1. for each file, get the full spect
% 2. calc the median and spread
% 3. go through each call and define the blobs that are large and exceed
% the threshold
% 4. if there are none, or they are too small, remove the call.

%filename='G:\USV data\Detections\C1_P6_1BL 2021-10-07  1_56 PM.mat';

%{
mydir='G:\USV data\Detections'; %uigetdir([],'Get folder with all the raw call mats');
allfiles=getAllFiles(mydir,'.mat');

wb=waitbar(0,'Starting up the curation process');
for fi=1:length(allfiles)

    load(allfiles{fi});
    % only go to 3 minutes
    % first only grab calls that start within 180 seconds
    Calls=Calls(Calls.Box(:,1)<180,:);
    
    if ~ismember('blobs',Calls.Properties.VariableNames)
        
        Calls.blobs{1}=struct([]);
        gpuDevice(1); % reset the gpu by clearing it
        audiodata.Data=gpuArray(audioread(audiodata.Filename,[1 min([180 audiodata.Duration])*audiodata.SampleRate])); % pull all samples
        
        % now load the full spectrogram onto gpu
        audiodata.windowsize = round(audiodata.SampleRate*0.0016); % hardcoded in for large spectrogram at least
        audiodata.noverlap = round(audiodata.SampleRate*.0008);
        audiodata.nfft = round(audiodata.SampleRate*.0016);
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
  

%}

%%
% try USVSEG filtering method


%load(filename);

%wav=audioread(audiodata.Filename,[1 min([180 audiodata.Duration])*audiodata.SampleRate]);
% first generate a spectrogram
%this takes wayyyyyy too long

%{
% the problem wit hthis
% initiate params struct
prm=struct('fs',audiodata.SampleRate);

% most of these parameters in table 2 in tachibanu 2020
prm.freqmin = 15000; prm.freqmax = 90000; % map freqs
prm.threshval = 4.5;   % std thresh for rat they suggest between 3 and 6 for rat
prm.durmin = 0.003;  prm.durmax = 1000; % for rat they suggest between 500 and 3k, min 3ms
prm.gapmin = 0.030; % min gap 30 to 40   
prm.margin = 0.015;  
prm.readsize = 10;

% now the multitaper spectrogram parameters
prm.fftsize=512; % datapoints in the fft (so quite a fat window)
prm.timestep=.0002; % 1/5 millisecond
prm.ntapers = 6; % 6 tapers
prm.step = round(prm.timestep*prm.fs); % this is the step single step (~70 datapts) 
prm.wavlen = length(wav); % should be around 70 million in my data (320k * 180)
prm.tapers = dpss(prm.fftsize,3,prm.tapers); % 6 tapers, time band prod is 3 (3 cycles in the taper)

prm.nsteps = floor((prm.wavlen-prm.fftsize+prm.step)/prm.step); % n steps will be done
prm.spgsize = prm.fftsize/2+1; % this is the frequencies i think
idx = repmat((1:prm.fftsize)',1,prm.nsteps)+repmat((0:prm.nsteps-1)*prm.step,prm.fftsize,1);
wavslice = wav(idx);
spmat = zeros(prm.spgsize,prm.nsteps);
parflag=1; prog=zeros(prm.nsteps,1);
if parflag == 1
    % this is insanely slow...
    parfor i=1:size(idx,2)
        spmat(i,:,:) = fft(wavslice(:,i).*prm.tapers,prm.fftsize); % fft all columns (each  of 6 tapers)

    end
else
    for n=1:ntapers
        spmat = fft(wavslice.*repmat(tapers(:,n),1,nsteps),fftsize);
    end
end
% average, and 
mtsp = 20*log10(squeeze(mean(spmat,2))*sqrt(1/(2*pi*fftsize)));
ft(1:(fftsize/2+1),:)
fvec = [0; (1:(prm.fftsize/2))'/prm.fftsize*prm.fs];
tvec = ((0:(size(mtsp,2)-1))*step+fftsize/2)/fs;

figure; imagesc(zscore(log(mean(abs(spmat(1:10000,1:(prm.fftsize/2+1),:)),3)))')
set(gca,'YDir','normal')
%}

%% This builds the spect with DeepSqueak code, and then cleans it using
% USVseg code.  USVseg is superior in cleaning spectrograms

% deepsqueak uses a quick spectrogram function that is based off of a fft
% usvseg writes their own multitaper function that basicaly is the above,
% but each short time fft also has many many tapers.



%{
% we can get a very high res
figure;
span=audiodata.SampleRate*10; % 10 seconds of data
sp=subplot(1,2,1);


% we can afford a really high res spectrogram (half the window, twice the
% freqs)
[ssx, fvec, tvec] = spectrogram(wav(1:span),round(1229/4),round(614/4),1229*2,audiodata.SampleRate,'yaxis');
spect=log(abs(ssx));

imagesc(tvec,fvec,spect);
set(gca,'YDir','normal','YLim',[15000 100000],'CLim',[-5 1])

% higher res is better
sp(2)=subplot(1,2,2);
[spect2, fvec, tvec] = spectrogram(wav(1:span),1229,614,1229,audiodata.SampleRate,'yaxis');
imagesc(tvec,fvec,log(abs(spect2)));
set(gca,'YDir','normal','YLim',[15000 100000])
linkaxes(sp);



% lifter
flims=[15000 90000];
figure; sp=subplot(1,3,1);
imagesc(tvec,fvec,spect);
set(gca,'YDir','normal','Clim',[-5 1],'ylim',flims);


% liftering the data
liftercutoff = 3; % fixed parameter
fftsize = (size(spect,1)-1)*2;
cep = fft([spect;flipud(spect(2:end-1,:))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);

sp(2)=subplot(1,3,2);
imagesc(tvec,fvec,log(abs(liftered)));
set(gca,'ydir','normal','Clim',[-1 3],'ylim',flims);
linkaxes(sp);
sp(3)=subplot(1,3,3);
imagesc(tvec,fvec,liftered);
set(gca,'ydir','normal','Clim',[-1 4],'ylim',flims);
linkaxes(sp);


% Median centering
med = median(liftered,2);
liftmed = liftered-repmat(med,1,size(liftered,2));
figure;
sp=subplot(1,2,1);
imagesc(tvec,fvec,liftered);
set(gca,'YDir','normal','Clim',[-1 3],'ylim',flims);

%set(gca,'ydir','normal','Clim',[-5 1]);
title('raw');
sp(2)=subplot(1,2,2);
imagesc(tvec,fvec,liftmed);
set(gca,'ydir','normal','Clim',[-1 3],'ylim',flims);
linkaxes(sp);
title('liftered median centered');

% movmedian smoothed

% i dont think we need this because were already at a high resolution
% we could highpass in the y axis though...


if exist('movmedian')==5
    fltnd = movmedian(liftmed,5,1);
else % for R2015b or earlier
    pad1 = repmat(liftmed(1,:),2,1);
    pad2 = repmat(liftmed(end,:),2,1);
    padded = [pad1;liftmed;pad2];
    fltnd = zeros(size(liftmed));
    for n=1:size(liftmed,1)
        fltnd(n,:) = median(padded(n+(0:4),:));
    end
end

fltnd=liftmed-movmean(liftmed,30,1);

figure;
sp=subplot(1,3,1);
imagesc(tvec,fvec,spect);
set(gca,'YDir','normal','Clim',[-5 1],'ylim',flims);
title('raw');

sp(2)=subplot(1,3,2);
imagesc(tvec,fvec,liftmed);
set(gca,'YDir','normal','Clim',[-1 3],'ylim',flims);
title('liftered median centered median smoothed');

sp(3)=subplot(1,3,3);
imagesc(tvec,fvec,fltnd);
set(gca,'YDir','normal','Clim',[-1 3],'ylim',flims);
title('thresholded');
linkaxes(sp);

% threshold


figure;
sp=subplot(1,2,1);
imagesc(tvec,fvec,spect);
set(gca,'YDir','normal','Clim',[-5 1],'ylim',flims);
title('raw');

sp(2)=subplot(1,2,2);
imagesc(tvec,fvec,zscore(liftmed,1,1)>2.5);
set(gca,'YDir','normal','ylim',flims);
title('liftered median centered median smoothed');
linkaxes(sp);

%}
%%
% here you can play with the fspan and the window/stepsize
% now toy with function
%{
[spect, f, t] = spectrogram(wav,1229,614,1229,audiodata.SampleRate,'yaxis');

mtsp=abs(spect).*sqrt(f);

[smoothed,centered]=cleanVoiceSpect(mtsp);

figure; 
sp=subplot(3,1,1);
imagesc(tvec(span),fvec,log(abs(mtsp)));
set(gca,'YDir','normal','CLim',[0 6]);
sp(2)=subplot(3,1,2);
imagesc(tvec(span),fvec,log(abs(centered)));
set(gca,'YDir','normal','Clim',[0 6]);
sp(3)=subplot(3,1,3);
imagesc(tvec(span),fvec,log(abs(smoothed)));
linkaxes(sp,'x','y');
set(gca,'YDir','normal','CLim',[-0 6]);
%}

