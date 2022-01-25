function [adjusted,fvec,tvec,liftmed,liftered,spect,ssx] = usvsegSpect(audiodata,params)
% This is the sequence USVseg uses to generate a flattened, adjusted
% spectrogram


try
    temp=audioinfo(audiodata.Filename);
catch
    [myfile,mydir]=uigetfile('',sprintf('Find audiofile for %s',audiodata.Filename));
    audiodata.Filename=fullfile(mydir,myfile);
end

if nargin<2
    params.winsize=.0016;
    params.fvec=0:120:100000;
end

winbins=round(audiodata.SampleRate*params.winsize);

stepsize=round(winbins/2); % use a 2x overlap

nfft=1229*2; % use double, you kill the top half anyways
% or use fvec....

gpuDevice(1);
% need to divide tbins into roughly 30 second bins
nsegments=ceil(diff(params.timewin)/30);
tbins=linspace(0,min([180 audiodata.Duration]),nsegments+1);
tbins=[tbins(1:end-1)' tbins(2:end)'];
tbins=round(tbins*audiodata.SampleRate); tbins(:,1)=tbins(:,1)+1;

useGPU=1;
if useGPU % need to figure out how to use the gpu on ths
    gpuDevice(1); % clear gpu
    spect=[]; tvec=[]; addon=0;% will need to preallocate this space later
    % parse into minute segments
    for i=1:6 % 3 minute recordings
        wav=gpuArray(audioread(audiodata.Filename,tbins(i,:)));
        [ssx, ftemp, ttemp] = spectrogram(wav,winbins,stepsize,params.fvec,audiodata.SampleRate,'yaxis');
        spect=[spect gather(10*log(ssx))];
        tvec=[tvec gather(ttemp)+addon(end)];
        addon=[addon max(tvec)];
    end
    fvec=gather(ftemp);
else
    wav=audioread(audiodata.Filename,[1 min([180 audiodata.Duration])*audiodata.SampleRate]);
    [ssx, fvec, tvec] = spectrogram(wav,winbins,stepsize,params.fvec,audiodata.SampleRate,'yaxis');
    spect=10*log(ssx);
end

% for some reason this creates oscillations in the spect image in the freq
% domain
adjusted = imadjust(imcomplement(abs(spect)./max(abs(spect(:)))));


% liftering the data- basically removes any slow oscillations in the
% frequency spectrum (like a standing oscillation, or any harmonics)
liftercutoff = 3; % fixed parameter
fftsize = (size(adjusted,1)-1)*2;
cep = fft([adjusted;flipud(adjusted(2:end-1,:))]);
lifter = ones(size(cep));
lifter(1:liftercutoff,:) = 0;
lifter((fftsize-liftercutoff+1):fftsize,:) = 0;
temp = real(ifft(cep.*lifter));
liftered = temp(1:(fftsize/2+1),:);



% median centering
med = median(liftered,2);
liftmed = liftered-repmat(med,1,size(liftered,2));

% smooth in y dir
% apriori i want to smooth over 1 msec and 1000 hz
%dimensions = [4 1];
%se = strel('rectangle', dimensions);
%dilated=imdilate(liftmed,se); % either dilate, or 
dilated=SmoothMat2(liftmed,[10 10],[1/1.6 1000/120]); 
%overdilated=filter([1 1 1 0 0 0 0 1 1 1]',1,dilated2);

%{
figure; 
sp=subplot(6,1,1); imagesc(abs(spect(:,1:100000)));
sp(2)=subplot(6,1,2);  imagesc(adjusted(:,1:100000));
sp(3)=subplot(6,1,3); imagesc(liftered(:,1:100000));
sp(4)=subplot(6,1,4);  imagesc(liftmed(:,1:100000));
sp(5)=subplot(6,1,5); imagesc(dilated(:,1:100000));

linkaxes(sp);
drawnow;
%}
adjusted=dilated;

% zscoring cleans up during call but messies between calls.
%adjusted=zscore(dilated,1,1); 




end

