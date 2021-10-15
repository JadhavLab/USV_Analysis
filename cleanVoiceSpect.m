function [smoothed,medcentered,liftered] = cleanVoiceSpect(spect,liftercutoff,movingspan)
% [smoothed,medcentered,liftered] = cleanVoiceSpect(spect,liftercutoff,movingspan)
% this function lifters, then median centers, then smooths a dataset
% The spectrum has to come in as the real component of a mtspectrogram or a
% fft spectrogram and it cant be logged.

if ~exist('liftercutoff','var') || (isempty(liftercutoff))
    liftercutoff = 3; % fixed parameter
end

% lifter
fftsize = (size(spect,1)-1)*2;
cep = fft([spect;flipud(spect(2:end-1,:))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);

% median filter
med = median(liftered,2);
medcentered = liftered-repmat(med,1,size(liftered,2));

if ~exist('movingspan','var') || (isempty(movingspan))
    movingspan = 5; % fixed parameter, in my case its ~.1 msec
end


% smooth with very short span
smoothed = movmedian(medcentered,movingspan);

% now subtract a highly smoothed image, like way high
oversmoothed=SmoothMat2(smoothed,[50 50],10);

%highpass=smoothed-oversmoothed;



end

%%
%{
figure; 
sp=subplot(3,1,1);
imagesc(smoothed); set(gca,'Ydir','normal','CLim',[0 100]);
sp(2)=subplot(3,1,2);
imagesc(oversmoothed); set(gca,'Ydir','normal','CLim',[0 5]);

sp(3)=subplot(3,1,3);
imagesc(smoothed-oversmoothed); set(gca,'Ydir','normal','CLim',[0 100]);
linkaxes(sp);
%}