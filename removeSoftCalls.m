%%removeSoftCalls



%%%% reject calls that are too soft %%%%%%%%%

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
%%
% try his filtering method
filename='G:\USV data\Detections\C1_P6_1BL 2021-10-07  1_56 PM.mat';
load(filename);

wav=audioread(audiodata.Filename,[1 min([180 audiodata.Duration])*audiodata.SampleRate]);
% first generate a spectrogram

fs=audiodata.SampleRate;
fftsize=512; % datapoints
timestep=.0002; % 1/5 millisecond
step = round(timestep*fs);
wavlen = length(wav);
tapers = dpsstapers;
ntapers = size(tapers,2);
nsteps = floor((wavlen-fftsize+step)/step);
spgsize = fftsize/2+1;
idx = repmat((1:fftsize)',1,nsteps)+repmat((0:nsteps-1)*step,fftsize,1);
wavslice = wav(idx);
spmat = zeros(spgsize,nsteps);
parflag=1;
if parflag == 1
    % use "parfor" if Parallel Computing Toolbox installed
    parfor n=1:ntapers 
        ft = fft(wavslice.*repmat(tapers(:,n),1,nsteps),fftsize);
        spmat = spmat + abs(ft(1:(fftsize/2+1),:));
    end
else
    for n=1:ntapers
        ft = fft(wavslice.*repmat(tapers(:,n),1,nsteps),fftsize);
        spmat = spmat + abs(ft(1:(fftsize/2+1),:));
    end
end
mtsp = 20*log10((spmat/ntapers)*sqrt(1/(2*pi*fftsize)));

[spect, f, t] = spectrogram(wav,1229,614,1229,audiodata.SampleRate,'yaxis');

mtsp=abs(spect).*sqrt(f);

fvec = [0; (1:(fftsize/2))'/fftsize*fs];
tvec = ((0:(size(mtsp,2)-1))*step+fftsize/2)/fs;

span=1:find(tvec<10,1,'last'); % this is 10 seconds of data
figure;
sp=subplot(2,1,1);
imagesc(tvec(span),fvec,log(mtsp(:,span)));

%set(gca,'ydir','normal','Clim',[-100 -50]);

%set(gca,'ydir','normal','Clim',[-1 5]);


% liftering the data
liftercutoff = 3; % fixed parameter
fftsize = (size(mtsp(:,span),1)-1)*2;
cep = fft([mtsp(:,span);flipud(mtsp(2:end-1,span))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);
sp(2)=subplot(2,1,2);
imagesc(tvec(span),fvec,log(abs(liftered(:,span))));
%set(gca,'ydir','normal','Clim',[-5 1]);
linkaxes(sp);
% liftering seems to do nothing demonstrable to the data

% Median centering


% median centering
med = median(liftered,2);
liftmed = liftered-repmat(med,1,size(liftered,2));
figure;
sp=subplot(2,1,1);
imagesc(tvec(span),fvec,mtsp(:,span));
%set(gca,'ydir','normal','Clim',[-5 1]);
title('raw');
sp(2)=subplot(2,1,2);
imagesc(tvec(span),fvec,liftmed(:,span));
%set(gca,'ydir','normal','Clim',[-5 1]);
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
sp=subplot(3,1,1);
imagesc(tvec(span),fvec,mtsp(:,span));
%set(gca,'ydir','normal','Clim',[-5 1]);
title('raw');
sp(2)=subplot(3,1,2);
imagesc(tvec(span),fvec,fltnd);
title('liftered median centered median smoothed');

sp(3)=subplot(3,1,3);
imagesc(tvec(span),fvec,fltnd>10);
%set(gca,'ydir','normal','Clim',[-5 1]);
linkaxes(sp);
title('thresholded');
%%

% now toy with function

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


