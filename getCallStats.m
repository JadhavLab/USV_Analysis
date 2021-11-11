function [Calls,oldCalls] = getCallStats(Calls,params,spect,thrshd)
% CallStats=getCallStats(Calls,spect,thrshd);
%  


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
%}

cstats=table([],'VariableNames',{'nSyll'});

histbins=15000:2000:90000;
freqhist=histbins(1:end-1)+mean(diff(histbins))/2;

durmin=round(.002/params.timestep);
tstep=mean(diff(params.tvec));
fstep=mean(diff(params.fvec));


% for each call:
for id=1:height(Calls)

    onsetI=max([1 find(params.tvec<=Calls.onsetTime(id,1),1,'last')]);
    offsetI=min([find(params.tvec>Calls.offsetTime(id,1),1,'first') length(params.tvec)]);
    callThrsh=logical(thrshd(:,onsetI:offsetI));
    callSpect=spect(:,onsetI:offsetI);
    % first gather the blobs from this image:
     callblobs=regionprops(callThrsh,callSpect,...
        'PixelIdxList','PixelList','Area','BoundingBox','MeanIntensity',...
        'Centroid','Orientation','Eccentricity','Extent','EulerNumber','FilledImage');
    % has to have at least
    % blobs have to be:
    % at least 2 msec, should be about 5 pix across

    % mean intensity must be above threshold of 2.2 sd above mean
    % orientation is greater than say 85*, or basically vertical
    % must not be 'holey' e.g. filled image cant be more than 110% of
    % holes image
    % most have a centroid higher than 180kHz
    % blob must not be through start or end of the call box

    okidx= cellfun(@(a) a(3)>durmin, {callblobs.BoundingBox}) &...
        [callblobs.MeanIntensity] > 2.2 & ...
        abs([callblobs.Orientation]) < 85 &...
        (cellfun(@(a) sum(a(:)), {callblobs.FilledImage})./[callblobs.Area]<1.1) & ...
        cellfun(@(a) a(1)>5 && (a(1)+a(3))<(offsetI-onsetI-5),{callblobs.BoundingBox});
    callblobs=callblobs(okidx);


     % reconstitute the blobs only:
    callBW=zeros(size(thrshd(:,onsetI:offsetI)));
    callIM=nan(size(thrshd(:,onsetI:offsetI)));
    for bl=1:length(callblobs)
        inBox=ceil(callblobs(bl).BoundingBox);
        callBW(inBox(2):inBox(2)+inBox(4)-1,inBox(1):inBox(1)+inBox(3)-1)=...
            callblobs(bl).FilledImage*bl;
        callIM(callblobs(bl).PixelIdxList)=callSpect(callblobs(bl).PixelIdxList);
    end

    allcallidx=cell2mat({callblobs.PixelList}');


    % now that we have our legit blobs, lets get some stats here...

    cstats.maxFreq(id)=params.fvec(find(sum(callBW,2)>0,1,'last')); % max frequency 
    cstats.minFreq(id)=params.fvec(find(sum(callBW,2)>0,1,'first')); % min frequency
    cstats.FreqSpread(id)=minFreq-maxFreq; % freq spread in hz
    [histy]=histcounts(params.fvec(allcallidx(:,2)),histbins); % spectrogram overall, normalized to 0-1
    cstats.specGram(id)=histy/max(histy); % this lets secondary peaks shine
    [~,modeind]=max(histy); 
    cstats.freqMode=freqhist(modeind); % histogram freq that is most common

    [cstats.maxIntensity(id),maxpos]=max(callIM(:)); % max intensity, time pos of max intensity
    cstats.maxIntensityF(id)=params.fvec(rem(maxpos,size(callIM,1))); % freq at max intensity
    cstats.maxIntensityT(id)=floor(maxpos/size(callIM,1))/size(callIM,2); % %% of time of call at max intensity
    cstats.MeanIntensity(id)=mean(callIM(:),'omitnan'); % mean intensity of spline

    % build a spline
    [~,splineinds]=max(callIM); 
    splinefreqs=params.fvec(splineinds);
    splinefreqs(splinefreqs==params.fvec(1))=nan;
    mdl=fitlm(1:sum(~isnan(splinefreqs))*tstep,splinefreqs(~isnan(splinefreqs)));
    
    callSlope=mdl.Coefficients.Estimate(2);
    callCenter=mdl.Coefficients.Estimate(1);

    % max instantaneous slope, max jump slope, mean freq dev (in spline)
    % mean abs freq dev in spline
    cstats.maxInstaSlope(id)=max(diff(splinefreqs)); % max insta freq spread (overall)
    cstats.maxJumpSlope(id)=max(diff(splinefreqs(~isnan(splinefreqs)))); % max jump slope(including jumps)
    cstats.firstMoment(id)=mean(diff(splinefreqs),'omitnan'); % mean freq delta
    cstats.firstAbsMoment(id)=mean(abs(diff(splinefreqs)),'omitnan'); % mean abs freq mod
    cstats.secondAbsMoment(id)=mean(abs(diff(diff(splinefreqs))),'omitnan'); % diff diff freq mod
    
    
    splineReal=splineinds; splineReal(splineReal==1)=nan; % mea
    [~,cstats.maxSplineFreq]=params.fvec(max(splineinds)); % max spline freq
    [~,cstats.minSplineFreq]=params.fvec(min(splineReal)); % min spline freq


    % now per blob
    % max blob mean slope, min blob mean slope, max abs slope, mean in
    % syllable spread (f)
    cstats.maxSyllSlope(id)=max(cell2mat({callblobs.Orientation})); % max within syll slope
    cstats.minSyllSlope(id)=min(cell2mat({callblobs.Orientation})); % min within syll slope (overall slope)
    cstats.maxAbsSyllSlope(id)=max(abs(cell2mat({callblobs.Orientation}))); % max abs blob orientation
    
    % create coords
    bitCoords=ceil(cell2mat({callblobs.BoundingBox}'));
    bitLengths=bitCoords(:,3)-1; 
    cstats.longestSyll(id)=max(bitLengths)/(offsetI-onsetI); % as a %% of time
    cstats.shortestSyll(id)=min(bitLengths)/(offsetI-onsetI); % as a %% of time
    cstats.longestSyllSec(id)=max(bitLengths)*tstep; % in seconds
    cstats.shortestSyllSec(id)=min(bitLengths)*tstep; % in seconds
    
    cstats.nSyll(id)=size(callblobs,1);
    cstats.ngaps(id)=sum(diff(splineinds>1)~=0)/2-1; % number of gaps inbetween call sylls
    
    % I see alot of harmonics that are just overlapping the lower freq, how
    % can I parse those?  The issue is that these may artifically suggest a
    % simple call is multiple...
    % maybe n nonoverapping syllables
    % first get overlapping areas
    % spline each wave:
    syllSpline=zeros(size(callIM,1),size(callIM,2));
    thisim=callIM; thisim(isnan(thisim))=0; myvec=1:size(callIM,2);
    mySyll=table([],[],[],[],'VariableNames',{'Slope','Var','Concavity','Curvature'});
    for sl=1:max(callBW(:))
        [~,myspline]=max((callBW==sl).*thisim);
        theseinds=sub2ind([size(thisim,1) size(thisim,2)],myspline(myspline>1),myvec(myspline>1));
        syllSpline(theseinds)=sl;
        % and also calculate the linear regression of these vals
        D=pdist([myvec(myspline>1)',myspline(myspline>1)'],'Euclidean');
        Z = squareform(D);
        leng=Z(1,end);
        c=0;
        for ll=2:length(Z)
            c=c+1;
            totleng(c)=Z(ll-1,ll);
        end
        cstats.Sinuosity(id)=sum(totleng)/leng;
        mySyllVar(sl)=var(params.fvec(myspline(myspline>1)),'omitnan'); % in Hz 
        
        mdl=fitlm(myvec(myspline>1)*tstep,myspline(myspline>1));
        mySyll.Slope(sl)=mdl.Coefficients.Estimate(1);
        mySyll.Concavity(sl)=myspline(myspline>1)-mdl.Fitted;
        mySyll.Curvature(sl)=mdl.MSE;
    end
    cstats.maxSyllVar(id)=max(mySyllVar); cstats.minSyllVar(id)=min(mySyllVar);
    cstats.MaxSyllSlope=max(mySyll.Slope);
    cstats.minSyllSlope=min(mySyll.Slope);
    
    cstats.minSyllSlope=max(abs(mySyll.Slope));
    cstats.maxSyllCurvature=max(mySyll.Curvature); % mse
    
    cstats.maxSyllConcavity=max(mySyll.Concavity); % signed
    cstats.minSyllConcavity=min(mySyll.Concavity); % signed
    
    cstats.meanSyllSlope=mean(mySyll.Slope);
    cstats.meanSyllConcavity=mean(mySyll.Concavity);

    % concavity (this is a stretch but i think i can do it
    % split call in half, see difference in mean slope part 1 vs part 2
    % now calculate curvature of this spline
    


    % %% of call that is actually voiced
    cstats.deadSpace(id)=sum(bitCoords(:,3)-1)/(bitCoords(end,1)+bitCoords(end,3)-1-bitCoords(1));
    % actual seconds of dead space
    cstats.deadSpaceSec(id)=((bitCoords(end,1) +bitCoords(end,3)-1-bitCoords(1))...
        -sum(bitCoords(:,3)-1))*tstep;
    if length(callblobs)>1
        syllCenters=cell2mat({callblobs.Centroid}');
        cstats.maxCrossSyllSpread(id)=max(max(abs(syllCenters(:,2)-syllCenters(:,2)')));
        cstats.crossSyllVar(id)=var(linearize(abs(syllCenters(:,2)-syllCenters(:,2)')),'omitnan');
    else
        cstats.maxCrossSyllSpread(id)=0;
        cstats.crossSyllVar(id)=0;
    end


    % tallest syllable

    
    % largest instantaneous freq step
    % steepest blob, shallowest blob, mean blob steepness
    % peak position in time (0-1 beginning to end)
    % valley position in time (0-1) % e.g. lowest freq in time
    % Freq distrib weighted by amplitude... just mean across y axis
    % something like concavity- like corr yvals from center out
    % need to remove harmonics, so find any blobs that have overlapping x,
    % and use the blob with the lowest slope (noise is vertical, harmonics
    % should be twice the slope)
    % spline derivitive (maybe include nans here), second deriv too(abs
    % tho) e.g. nonlinear energy
    % area under the curve, so take the lowest freq, make that zero, then
    % pull the area under the curve of that. do this in absolute hz not
    % relative, and could get area fo reach blob
    




%
end