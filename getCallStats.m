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




% for each call:
for id=1:height(Calls)

    onsetI=max([1 find(params.tvec<=Calls.Box(id,1),1,'last')]);
    offsetI=min([find(params.tvec>Calls.Box(id,1)+Calls.Box(id,3),1,'first') length(params.tvec)]);
    callThrsh=logical(thrshd(:,onsetI:offsetI));
    callSpect=spect(:,onsetI:offsetI);
    % first gather the blobs from this image:
     callblobs=regionprops(callThrsh,callSpect,...
        'PixelIdxList','PixelList','Area','BoundingBox','MeanIntensity',...
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
    callBW=zeros(size(thrshd(:,onsetI:offsetI)));
    callIM=nan(size(thrshd(:,onsetI:offsetI)));
    for bl=1:length(callblobs)
        inBox=ceil(callblobs(bl).BoundingBox);
        callBW(inBox(2):inBox(2)+inBox(4)-1,inBox(1):inBox(1)+inBox(3)-1)=...
            callblobs(bl).FilledImage;
        callIM(callblobs(bl).PixelIdxList)=callSpect(callblobs(bl).PixelIdxList);
    end

    allcallidx=cell2mat({callblobs.PixelList}');


    % now that we have our legit blobs, lets get some stats here...
    maxFreq=params.fvec(find(sum(callBW,2)>0,1,'last'));
    minFreq=params.fvec(find(sum(callBW,2)>0,1,'first'));
    [histy]=histcounts(params.fvec(allcallidx(:,2)),15000:2000:90000); % spectrogram overall, normalized to 0-1
    specGram=histy/max(histy); % this lets secondary peaks shine

    [maxIntensity,maxpos]=max(callIM(:));
    maxIntensityF=params.fvec(rem(maxpos,size(callIM,1)));
    MeanIntensity=mean(callIM(:),'omitnan');

    % build a spline
    [~,splineinds]=max(callIM);
    splinefreqs=params.fvec(splineinds);
    splinefreqs(splinefreqs==params.fvec(1))=nan;

    % max instantaneous slope
    instaslope=max(diff(splinefreqs));
    jumpslope=max(diff(splinefreqs(~isnan(splinefreqs))));
    firstMoment=nanmean(diff(splinefreqs));
    firstAbsMoment=nanmean(abs(diff(splinefreqs)));


    % now per blob
    maxSlope=max(cell2mat({callblobs.Orientation}));
    minSlope=min(cell2mat({callblobs.Orientation}));
    
    maxAbsSlope=max(abs(cell2mat({callblobs.Orientation})));

    % create coords
    bitCoords=ceil(cell2mat({callblobs.BoundingBox}'));
    bitLengths=bitCoords(:,3)-1; 
    longestSyll=max(bitLengths)/(offsetI-onsetI); % as a %% of time
    shortestSyll=min(bitLengths)/(offsetI-onsetI); % as a %% of time
    nSyll=size(callblobs,1);


    % %% of call that is actually voiced
    deadSpace=sum(bitCoords(:,3)-1)/(bitCoords(end,1)+bitCoords(end,3)-1-bitCoords(1));

    syllCenters=cell2mat({callblobs.Centroid}');
    maxSyllSpread=max(max(abs(syllCenters(:,2)-syllCenters(:,2)')));

    SyllVar=var(linearize(abs(syllCenters(:,2)-syllCenters(:,2)')),'omitnan');

    % now need to do sinuosity


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