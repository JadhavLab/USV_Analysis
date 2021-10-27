%usvsegDetect

%% Make the spectrogram


load('G:\USV data\Detections\C3-P20-T8 2021-10-07  3_26 PM.mat');
% variables are audiodata, Calls, detection_metadata

winsize=.0008;
winbins=round(audiodata.SampleRate*winsize);

stepsize=round(winbins/2); % use a 2x overlap

nfft=1229*2; % use double, you kill the top half anyways
% or use fvec....
fvec=0:120:100000;

% spectrogram settings
params=struct('fs',audiodata.SampleRate);
params.winsize = 0.0008;  params.fvec=0:120:100000;
params.timestep=params.winsize/2; % will always round when converting to indices

% call settings
params.freqmin = 15000; params.freqmax = 120000;
params.threshval = 2.1;   params.durmin = 0.002;  params.durmax = 0.6;
params.gapmin = 0.030;    params.margin = params.gapmin/2;  

% this bit takes FOREVER, because im using very hgih res 
[adjusted,fvec,tvec] = usvsegSpect(audiodata,params);

% process all the way- get the segments here
[spect,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq,thresh,contflg] =...
    usvsegProcfun(adjusted,params);

% audiodata is how to get the data...
% so we need to whiten the whole dataset.  to do this, we zscore each freq
% independently.  Luckily i have a gpu that can handle that...


threshval=2.45; % sd above the mean (this could be 2) % or 2/100


%%
