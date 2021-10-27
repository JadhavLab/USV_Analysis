function [adjusted,fvec,tvec,liftmed,liftered,spect,ssx] = usvsegSpect(audiodata,params)
% This is the sequence USVseg uses to generate a flattened, adjusted
% spectrogram


if nargin<2
    params.winsize=.0008;
    params.fvec=0:120:100000;
end

winbins=round(audiodata.SampleRate*params.winsize);

stepsize=round(winbins/2); % use a 2x overlap

nfft=1229*2; % use double, you kill the top half anyways
% or use fvec....


wav=audioread(audiodata.Filename,[1 min([180 audiodata.Duration])*audiodata.SampleRate]);

[ssx, fvec, tvec] = spectrogram(wav,winbins,stepsize,params.fvec,audiodata.SampleRate,'yaxis');

spect=log(abs(ssx));



% liftering the data
liftercutoff = 3; % fixed parameter
fftsize = (size(spect,1)-1)*2;
cep = fft([spect;flipud(spect(2:end-1,:))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);



% median centering
med = median(liftered,2);
liftmed = liftered-repmat(med,1,size(liftered,2));


adjusted=zscore(liftmed,1,1);


end

