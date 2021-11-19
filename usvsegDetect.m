function [spect,thrshd,params,onoffset,onoffsetm,blobthrshd]=usvsegDetect(audiodata,timewin)
% function [spect,thrshd,onoffset,onoffsetm]=usvsegDetect(audiodata)
% function usvsegDetect turns the usvseg workflow into a quick function

% INPUTS:
%   audiodata: struct with necessary fields: SampleRate, Filename, and
%   Duration
% OUTPUTS:
%   Spect: a corrected spectrogram of the file
%   thrshd: a thresholded boolean matrix of the calls
%   onoffset: onoffsets of the ACTUAL CALLS
%   onoffsetm: onoffset with a margin so you can get before and after THE
%   call.

% Algorithm:
% 1.create a very high res fft spectrogram
%   .8 msec windowsize, half that for the bin (0.4msec) timestep
%   also high spectral resolution, about 120 hz per step.  so its already
%   very smooth in the frequency domain
% 2. Threshold
%   This uses a fwhm z-scoring method to capture variance without
%   considering the tail of the data
%   This also generates a corrected boolean 1's and zeros' matrix of
%   'loudness' calls
% 3. Find the ridges
%   Ridges have to be sufficiently large in Frequency and Time domains.
%   First is the frequency domain (opportunity to also use a saliency
%   metric here.
%   Second is the time domain- tries to comb over small gaps in the calls









%{
% test a file out here!
load('G:\USV data\Detections\C3-P20-T8 2021-10-07  3_26 PM.mat');
% variables are audiodata, Calls, detection_metadata
audiodata.Filename='E:\Brandeis datasets\FMR1 Project Data\USV data\raw\Wav files\C3-P20-T8.wav';
%}
% initiate struct with parameters from audiofile


if ~exist('timewin','var') timewin=180; end % 3 minutes
if isempty(timewin), timewin=180; end
if length(timewin)==1, timewin=[0 timewin]; end % start to end
params=struct('fs',audiodata.SampleRate,'timewin',timewin);


% HARDCODED PARAMETERS:
winsize=.0008;
% bins in a window
winbins=round(audiodata.SampleRate*winsize);
% bins in a step
stepbins=round(winbins/2); % use a 2x overlap



% spectrogram settings
params.winsize = 0.0008;  params.fvec=0:300:100000;
params.timestep=params.winsize/2; % will always round when converting to indices

% call settings
params.freqmin = 15000; params.freqmax = 90000;
params.threshval = 2.2;   params.durmin = 0.003;  params.durmax = 0.6;
params.gapmin = 0.025;    params.margin = params.gapmin/2;  


% Make the spectrogram

% this bit takes FOREVER, because im using very hgih res 
[adjusted,fvec,tvec] = usvsegSpect(audiodata,params);

% add time and freq scales
params.fvec=fvec;
params.tvec=tvec;

% process all the way- get the segments here
[spect,thrshd,onoffset,onoffsetm,blobthrshd] =usvsegProcfun(adjusted,params);



% audiodata is how to get the data...
% so we need to whiten the whole dataset.  to do this, we zscore each freq
% independently.  Luckily i have a gpu that can handle that...
end



%%


% just to see what the spect looks like
%{

myinds=10000:20000;
figure; imagesc(tvec(myinds),fvec,adjusted(:,myinds));
% bracket y 
set(gca,'Ydir','normal','CLim',[-1 4],'YLim',[20000,100000]);
%}