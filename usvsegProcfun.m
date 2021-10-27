function [spect,onoffset,onoffsetm,freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq,thresh,contflg] = usvsegProcfun(spect,params)

fs=params.fs;
timestep = params.timestep;
gapmin = params.gapmin;
durmin = params.durmin;
durmax = params.durmax;
margin = params.margin;
freqmin = params.freqmin;
freqmax = params.freqmax;


% threshold calculation with n*sigma (SD) of background noise 

thresh = estimatethresh(spect,fs,freqmin,freqmax,params.threshval);

% thresholding, thrshd is the thresholded image
thrshd = thresholding(spect,fs,freqmin,freqmax,thresh);
% onset/offset detection
[onoffset,onoffsetm,~,contflg] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin);
% peak tracking
[freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq] = specpeaktracking(spect,fs,timestep,freqmin,freqmax,onoffset,margin);


end




% ////////////////////////////////////////////////////////////////////////
% this estimates the z-threshold by using the fwhm method of the
% distribution.  This method is better than the tail method because it will
% ignore outliers
function thresh = estimatethresh(fltnd,fs,freqmin,freqmax,threshval)
fftsize = (size(fltnd,1)-1)*2;
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
cut = fltnd(fnmin:fnmax,:);
bin = -0.05:0.1:10;
bc = bin(1:end-1)+diff(bin)/2;
h = histcounts(cut(:),bin);
fwhm = bc(find(h<h(1)/2,1)-1)*2;
sigma = fwhm/2.35;
thresh = sigma*threshval;

end


% ////////////////////////////////////////////////////////////////////////
% this threhsolds the sounds and also removes anythign that isnt in the
% search frequency by using a mask
function thrshd = thresholding(fltnd,fs,freqmin,freqmax,thresh)
fftsize = (size(fltnd,1)-1)*2;
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
mask = zeros(size(fltnd));
mask(fnmin:fnmax,:) = 1;
thrshd = (fltnd>thresh).*mask;
end







% ////////////////////////////////////////////////////////////////////////
% this detects
function [onoffset,onoffsetm,onoffsig,contflg] = detectonoffset(thrshd,fs,timestep,gapmin,durmin,durmax,margin,onoffthresh)
if nargin<8
    onoffthresh = 5; % optimized for multitaper spectrogram
end
fftsize = (size(thrshd,1)-1)*2;
step = round(timestep*fs);
% onset/offset detection
onoff = max(filter(ones(onoffthresh,1),1,thrshd))'>=onoffthresh;


% merge fragmented pieces
% this doesnt appear to do that.  It appears to shorten the duration of
% each and every call...
% 
ndurmin = round(durmin*fs/step);
f = filter(ones(ndurmin,1)/ndurmin,1,[onoff;zeros(round(ndurmin/2),1)]);
% you can add gap filters in...

monoff = f(round(ndurmin/2)+1:end)>0.5;
monoff(1) = 0; 
monoff(end) = onoff(end);
onidx = find(diff(monoff)>0)+1; %
offidx = find(diff(monoff)<0)+1; % 
if isempty(onidx)||isempty(offidx)
    onoffset = zeros(0,2);
    onoffsetm = zeros(0,2);
    onoffsig = zeros(size(thrshd,2),1);
    contflg = 0;
    return;
end
offidx(end) = min(offidx(end),length(monoff));
if ~isempty(onidx) && monoff(end)==1
    onidx = onidx(1:end-1);
end



% continuity flag: check if the end of read file is "ON"
contflg = monoff(end);

% gap thresholding
% if any gap between on-off pairs is too short, combine
gap = (onidx(2:end)-offidx(1:end-1))*timestep;
gap(end+1) = 0;
gid = find(gap>=gapmin);
if ~isempty(gid)
    onidx = [onidx(1); onidx(gid+1)];
    offidx = [offidx(gid); offidx(end)];
else
    onidx = onidx(1);
    offidx = offidx(end);
end

% syllable duration threholding
dur = (offidx-onidx)/fs*step;
did = find(durmin<=dur & dur<=durmax);
onidx = onidx(did);
offidx = offidx(did);
tvec = ((0:(size(thrshd,2)-1))'*step+fftsize/2)/fs;
onset = tvec(onidx);
offset = tvec(offidx);
if isempty(onset)||isempty(offset)
    onoffset = zeros(0,2);
    onoffsetm = zeros(0,2);
    onoffsig = zeros(size(thrshd,2),1);
    contflg = 0;
    return;
end
% margin addition
onsetm = onset-margin;
offsetm = offset+margin;
% syllables whose margins are overlapped are integrated in onoffsetm but not in onoffset
idx = find((onsetm(2:end)-offsetm(1:end-1))>0);
onsetI = [onset(1);onset(idx+1)];
offsetI = [offset(idx);offset(end)];
onsetm = onsetI-margin;
onsetm(1) = max(1/fs*step,onsetm(1));
offsetm = offsetI+margin;
offsetm(end) = min(max(size(thrshd,2)*step/fs),offsetm(end));    
% output 
onoffset = [onset offset];
onoffsetm = [onsetm offsetm];
% on/off signal
temp = zeros(size(onoff));
onidx2 = round((onset*fs-fftsize/2)/step+1);
offidx2 = round((offset*fs-fftsize/2)/step+1);
temp(onidx2) = 1;
temp(offidx2+1) = -1;
onoffsig = cumsum(temp); % throwaway variable, it has overlaps...
end


% ////////////////////////////////////////////////////////////////////////
function [freqtrace,amptrace,maxampval,maxampidx,maxfreq,meanfreq,cvfreq] = specpeaktracking(fltnd,fs,timestep,freqmin,freqmax,onoffset,margin)
% INPUTS
% fltnd= filtered image
% fs- samp freq
% timestep= timesteps
% onoffset=

if isempty(onoffset)
    freqtrace = nan(size(fltnd,2),4);
    amptrace =  nan(size(fltnd,2),4);
    maxampval = [];
    maxampidx = [];
    maxfreq = [];
    meanfreq = [];
    cvfreq = [];
    return;
end
fftsize = (size(fltnd,1)-1)*2;
step = round(fs*timestep);
% spectral peak saliency, basically filter to the frequency width we
% expect, you can either build this from measuring it, or guess it (i
% guessed it to be a sine wave about 50 wide, so you can filt with a sine
% and a cos to get the correct offset.
spcsal = spectralsaliency(fltnd);
% get peak and reorganize
bandwidth = 9;
ncandidates = 4;
contmin = 10;
ampthr = 0;

nstep = size(fltnd,2);
fnmin = floor(freqmin/fs*fftsize)+1;
fnmax = ceil(freqmax/fs*fftsize)+1;
onidx = round((onoffset(:,1)-margin)*fs/step); % add margin
offidx = round((onoffset(:,2)+margin)*fs/step); % add margin
onidx(1) = max(1,onidx(1));
offidx(end) = min(nstep,offidx(end));
freqmat = nan(nstep,ncandidates);
ampmat = nan(nstep,ncandidates);
maxampval = nan(length(onoffset(:,1)),1);
maxampidx = ones(length(onoffset(:,1)),1);
maxampfreq = nan(size(onoffset,1),1);
meanfreq = nan(size(onoffset,1),1);
cvfreq = nan(size(onoffset,1),1);

for n=1:size(onoffset,1)
    idx = onidx(n):offidx(n);
    [peakfreq,peakamp] = searchpeak(spcsal(fnmin:fnmax,idx),fltnd(fnmin:fnmax,idx),ncandidates,bandwidth);
    [peakfreqsg,peakampsg] = segregatepeak(peakfreq+fnmin-1,peakamp,contmin,ampthr);
    freqmat(idx,:) = peakfreqsg;
    ampmat(idx,:) = peakampsg;
    if any(~isnan(peakampsg(:)))
        [mvC,miC] = max(peakampsg,[],2);
        [~,miR] = max(mvC);
        maxampidx(n) = miR+idx(1)-1;
        maxampfreq(n) = peakfreqsg(miR,miC(miR));
        if ~isnan(maxampfreq(n))
            maxampval(n) = fltnd(round(maxampfreq(n)),maxampidx(n));
        end
    end
    meanfreq(n) = mean((freqmat(idx,1)-1)/fftsize*fs,'omitnan');
    ft = (peakfreqsg(:,1)-1)/fftsize*fs;
    cvfreq(n) = std(ft,1,'omitnan')/meanfreq(n);
end
freqtrace = (freqmat-1)/fftsize*fs;
amptrace = ampmat;
maxfreq = (maxampfreq-1)/fftsize*fs;
end


% ////////////////////////////////////////////////////////////////////////
function [peakfreq,peakamp] = searchpeak(specsaliency,specamp,ncandidates,bandwidth)
num_steps = size(specsaliency,2);
search_range = bandwidth-1;
remove_range = bandwidth*2-1;
peakfreq = nan(num_steps,ncandidates);
peakamp = nan(num_steps,ncandidates);
specsaliency(specsaliency<0) = 0;
for n=1:num_steps
    slice = specsaliency(:,n);
    for c=1:ncandidates
        [~,mi] = max(slice);
        % center of gravity
        rng = max(mi-search_range,1):min(mi+search_range,size(slice,1));
        temp = specsaliency(rng,n);
        peak = sum(temp.*rng')/sum(temp); 
        % store
        peakfreq(n,c) = peak;
        peakamp(n,c) = specamp(mi,n);
        % remove
        idx = max(round(peak)-remove_range,1):min(round(peak)+remove_range,size(slice,1));
        slice(idx) = -Inf;
    end
end

end
% ////////////////////////////////////////////////////////////////////////
function [peakfreqsg,peakampsg] = segregatepeak(peakfreq,peakamp,conthr,ampthr)
% amplitude thresholding
peakfreq(peakamp<ampthr) = nan;
peakamp(peakamp<ampthr) = nan;
%object segregatin with putting object number
%allow skipping two frames (steps)
distthr = 0.05; % 5 percent: fixed parameter
[nstep,ncand] = size(peakfreq);
objmat = reshape((1:(nstep*ncand)),ncand,nstep)';
objmat(isnan(peakfreq)) = nan;
nskip = 2; % can skip max 2 frames if intermediate framse are NaN
distmat = nan(nstep-3,ncand,nskip+1);
pathmat = nan(nstep-3,ncand,nskip+1);
for n=1:nstep-nskip-1
    for m=1:nskip+1
        temp = abs((1./peakfreq(n,:)'*peakfreq(n+m,:))-1);
        [mv,mid] = min(temp,[],2);
        distmat(n,:,m) = mv;
        pathmat(n,:,m) = mid;
    end
    pm = pathmat;
    pm(distmat>distthr) = nan;
    pm(isnan(distmat)) = nan;
    pp = pm(:,:,1);
    if any(~isnan(pm(n,:,1)))
        pp = pm(:,:,1);
        x = n+1;
    elseif any(~isnan(pm(n,:,2)))
        pp = pm(:,:,2);
        x = n+2;
    elseif any(~isnan(pm(n,:,3)))
        pp = pm(:,:,3);
        x = n+3;
    end
    for m=1:ncand
        if ~isnan(pp(n,m))
            if objmat(x,pp(n,m)) < objmat(n,m)
                val = objmat(x,pp(n,m));
                objmat(objmat==val) = objmat(n,m);
            else
                objmat(x,pp(n,m)) = objmat(n,m);
            end
        end
    end
end
% thresholding
objnum = unique(objmat(:));
objnum = objnum(~isnan(objnum));
peaks2 = peakfreq;
ampmat2 = peakamp;
objmat2 = objmat;
objlen = zeros(length(objnum),1);
objamp = zeros(length(objnum),1);
for n=1:length(objnum)
    idx = find(objmat==objnum(n));
    objlen(n) = length(idx);
    objamp(n) = mean(ampmat2(objmat==objnum(n)));
end
for n=1:length(objlen)
    if objlen(n)<conthr
        objlen(n) = nan;
        peaks2(objmat==objnum(n)) = nan;
        objmat2(objmat==objnum(n)) = nan;
        ampmat2(objmat==objnum(n)) = nan;
    end
end
objnum = objnum(~isnan(objlen));
objamp = objamp(~isnan(objlen));
objlen = objlen(~isnan(objlen));
% sorting
peakfreqsg = nan(size(peaks2));
peakampsg = nan(size(peakamp));
for n=1:nstep
    on = objmat2(n,:);
    oa = nan(length(on),1);
    for m=1:length(on)
        if ~isempty(find(objnum==on(m)))
            oa(m) = objamp(objnum==on(m));
        end
    end
    oa2 = oa;
    oa2(isnan(oa)) = -Inf;
    [~,sid] = sort(oa2,'descend');
    peakfreqsg(n,:) = peaks2(n,sid);
    peakampsg(n,:) = ampmat2(n,sid);
end

end

% ////////////////////////////////////////////////////////////////////////

function spcsal = spectralsaliency(fltnd)
% looks like it convolves with a wavelet
fftsize = (size(fltnd,1)-1)*2;
tfsp = fftshift(sum(abs(fft(dpsstapers,fftsize)),2));
dtfsp = -diff(tfsp,2); % second-order differential
rng = fftsize/2+(-6:6); % pull the center from the filters?
rdtfsp = dtfsp(rng); % now use the center of that centfreq function
%salfilt = (rdtfsp-mean(rdtfsp))/std(rdtfsp); % now convolve w that filter
% prob need to convolve with a lower freq filter

% i have a better filter...
salfilt=-cos(linspace(0,2*pi,50));
%salfilt=sin(linspace(0,2*pi,50));
filtsz=round(length(salfilt)/2);
fil = filter(salfilt,1,[fltnd;zeros(filtsz,size(fltnd,2))]);
spcsal = fil(filtsz+1:end,:);
end



