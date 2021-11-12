function [cstats,Calls] = getCallStats(Calls,params,spect,thrshd)
% CallStats=getCallStats(Calls,spect,thrshd);
%  


%{
First, lets understand the unput data that DeepSqueak uses
1. Entropy: geomean/mean of the full window thats just loudness
2. ridgeTime: the indices in the call box of the ridge (calced by either
	bs power or entropy (I use blobs instead)
3. ridgeFreq: they use a single frequency, i may be able to get away with a
    few
4. smoothed ridge frequency (using a rlowess smoothing algo, using .1 input
5. a filtered image using imgradientxy on the raw image I
6. signal to noise: mean 1-entropy of the ridge
7. begin time (of the ridge)
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

cstats=table(nan(height(Calls),1),'VariableNames',{'nSyll'});

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
    
    if isempty(callblobs)
        Calls.Accept(id)=0;
        %keyboard
        continue
    end

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
    cstats.nSyll(id)=size(callblobs,1);

    cstats.maxFreq(id)=params.fvec(find(sum(callBW,2)>0,1,'last')); % max frequency 
    [~,maxind]=max(find(callBW>0,1,'last')); % find last element in each col thats  1
    cstats.maxFreqInd=maxind/size(callBW,2);
    cstats.minFreq(id)=params.fvec(find(sum(callBW,2)>0,1,'first')); % min frequency
    [~,minind]=max(find(callBW>0,1,'first')); % find last element in each col thats  1
    cstats.minFreqInd=minind/size(callBW,2);
    cstats.FreqSpread(id)=cstats.maxFreq(id)-cstats.minFreq(id); % freq spread in hz
    [histy]=histcounts(params.fvec(allcallidx(:,2)),histbins); % spectrogram overall, normalized to 0-1
    cstats.specGram(id,:)=histy/max(histy); % this lets secondary peaks shine
    [~,modeind]=max(histy); 
    cstats.freqMode(id)=freqhist(modeind); % histogram freq that is most common

    [cstats.maxIntensity(id),maxpos]=max(callIM(:)); % max intensity, time pos of max intensity
    cstats.maxIntensityF(id)=params.fvec(rem(maxpos,size(callIM,1))); % freq at max intensity
    cstats.maxIntensityT(id)=floor(maxpos/size(callIM,1))/size(callIM,2); % %% of time of call at max intensity
    cstats.MeanIntensity(id)=mean(callIM(:),'omitnan'); % mean intensity of spline

    % build a spline
    [~,splineinds]=max(callIM); 
    splinefreqs=params.fvec(splineinds);
    splinefreqs(splinefreqs==params.fvec(1))=nan;
    mdl=fitlm((1:sum(~isnan(splinefreqs)))*tstep,splinefreqs(~isnan(splinefreqs))');
    
    cstats.callSlope(id)=mdl.Coefficients.Estimate(2);
    cstats.callCenter(id)=mdl.Coefficients.Estimate(1);

    % max instantaneous slope, max jump slope, mean freq dev (in spline)
    % mean abs freq dev in spline
    cstats.maxInstaSlope(id)=max(diff(splinefreqs)); % max insta freq spread (overall)
    cstats.maxJumpSlope(id)=max(diff(splinefreqs(~isnan(splinefreqs)))); % max jump slope(including jumps)
    cstats.firstMoment(id)=mean(diff(splinefreqs),'omitnan'); % mean freq delta
    cstats.firstAbsMoment(id)=mean(abs(diff(splinefreqs)),'omitnan'); % mean abs freq mod
    cstats.secondAbsMoment(id)=mean(abs(diff(diff(splinefreqs))),'omitnan'); % diff diff freq mod
    
    
    splineReal=splineinds; splineReal(splineReal==1)=nan; % mea

    cstats.maxSplineFreq(id)=params.fvec(max(splineinds));
    cstats.minSplineFreq(id)=params.fvec(min(splineReal));


    % now per blob
    % max blob mean slope, min blob mean slope, max abs slope, mean in
    % syllable spread (f)
    % lets use the real slope here not the oval
    %cstats.maxSyllSlope(id)=max(cell2mat({callblobs.Orientation})); % max within syll slope
    %cstats.minSyllSlope(id)=min(cell2mat({callblobs.Orientation})); % min within syll slope (overall slope)
    %cstats.maxAbsSyllSlope(id)=max(abs(cell2mat({callblobs.Orientation}))); % max abs blob orientation
    
    % create coords
    bitCoords=ceil(cell2mat({callblobs.BoundingBox}'));
    bitLengths=bitCoords(:,3)-1; 
    cstats.longestSyll(id)=max(bitLengths)/(offsetI-onsetI); % as a %% of time
    cstats.shortestSyll(id)=min(bitLengths)/(offsetI-onsetI); % as a %% of time
    cstats.longestSyllSec(id)=max(bitLengths)*tstep; % in seconds
    cstats.shortestSyllSec(id)=min(bitLengths)*tstep; % in seconds
    
    
    cstats.ngaps(id)=sum(diff(splineinds>1)~=0)/2-1; % number of gaps inbetween call sylls
    
    % I see alot of harmonics that are just overlapping the lower freq, how
    % can I parse those?  The issue is that these may artifically suggest a
    % simple call is multiple...
    % maybe n nonoverapping syllables
    % first get overlapping areas
    % spline each wave:
    syllSpline=zeros(size(callIM,1),size(callIM,2));
    thisim=callIM; thisim(isnan(thisim))=0; myvec=1:size(callIM,2);
    efill=nan(max(callBW(:)),1);
    mySyll=table(efill,efill,efill,efill,efill,...
        'VariableNames',{'Slope','Var','Concavity','Curvature','Height'});
    for sl=1:max(callBW(:))
        [~,myspline]=max((callBW==sl).*thisim);
        theseinds=sub2ind([size(thisim,1) size(thisim,2)],myspline(myspline>1),myvec(myspline>1));
        syllSpline(theseinds)=sl;
        % and also calculate the linear regression of these vals
        D=pdist([myvec(myspline>1)',myspline(myspline>1)'],'Euclidean');
        Z = squareform(D);
        leng=Z(1,end);
        c=0;
        totleng=[];
        for ll=2:length(Z)
            c=c+1;
            totleng(c)=Z(ll-1,ll);
        end
        cstats.Sinuosity(id)=sum(totleng)/leng;
        mySyll.Var(sl)=var(params.fvec(myspline(myspline>1)),'omitnan'); % in Hz 
        
        mdl=fitlm(myvec(myspline>1)*tstep,myspline(myspline>1));
        mySyll.Slope(sl)=mdl.Coefficients.Estimate(1);

        realspline=myspline(myspline>1);
        % inside quarters vs outside quarters
        quarters=round(linspace(1,4,length(realspline)));
        mySyll.Concavity(sl)=(sum(realspline(quarters==1 |quarters==4)'-mdl.Fitted(quarters==1 |quarters==4))-...
            sum(realspline(quarters==2 |quarters==3)'-mdl.Fitted(quarters==2 |quarters==3)))*fstep;
        mySyll.Curvature(sl)=mdl.MSE;
        mySyll.Height(sl)=max(fstep*(realspline-realspline'),[],'all');
    end
    cstats.maxSyllVar(id)=max(mySyll.Var); 
    cstats.minSyllVar(id)=min(mySyll.Var);

    cstats.maxSyllSlope(id)=max(mySyll.Slope);
    cstats.minSyllSlope(id)=min(mySyll.Slope);
    cstats.meanSyllSlope(id)=mean(mySyll.Slope);

    cstats.maxAbsSyllSlope(id)=max(abs(mySyll.Slope));
    cstats.maxSyllCurvature(id)=max(mySyll.Curvature); % mse
    
    cstats.maxSyllConcavity(id)=max(mySyll.Concavity); % signed
    cstats.minSyllConcavity(id)=min(mySyll.Concavity); % signed
    
    cstats.meanSyllSlope(id)=mean(mySyll.Slope);
    cstats.meanSyllConcavity(id)=mean(mySyll.Concavity);

    cstats.tallestSyll(id)=max(mySyll.Height);
    cstats.shortestSyll(id)=min(mySyll.Height);
    

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


%
end

cstats=cstats(Calls.Accept==1,:); Calls=Calls(Calls.Accept==1,:); 
%{
callstats2=cstats;
callstats2.specGram=[];
normalized=rescale(table2array(callstats2),'InputMin',...
    min(table2array(callstats2)),'InputMax',max(table2array(callstats2)));
figure; imagesc(normalized)
% lets see if these actually work okay

tree = linkage(normalized,'centroid');
figure()
dendrogram(tree,100)
% now lets see if we can cluster a single session a bit

for i=1:size(normalized,2)
    for j=i+1:size(normalized,2)-1
        figure;
        scatter(normalized(:,i),normalized(:,j))
        title(sprintf('%s %s',callstats2.Properties.VariableNames{i},callstats2.Properties.VariableNames{j}));
    end
end
%}



end