%{
This is the analysis sketchpad of the usv data from start to finish
 I will start with blocks of sketch code here and then offload them into
 functions or other scripts if they become useful.

the preprocessing pipeline was as follows:
-data were acquired, see acquisition protocol
-a systematically random selection of files were pulled and manually
curated for training data.
-a network was built on those calls using DeepSqueak v 3.0.1
-the full detection algorithm was then performed on the full gambit of data
using both the rat detector yolo r1 and rat detector JHB networks together
-a post-hoc denoiser was then applied to ALL FILES

Manual edits:
First some observations:

-the bounding boxes tend to miss the start and end of the calls, leaving
loud bits outside the box

- The contour isnt quite right sometimes, need to check how they threshold
to get a contour.
    The contour math:
    -they use a conjunctive approach.  they need at least 6 x coords in
    which... both the entropy (geomean/mean)>a threshold) and the loudness
    threshold are used. Loudness threshold=.825 and entropy is .215.. then
    they increase each by 10% until you get 6 dots
    -the elbow method could be used to find the top bit of spiking
    -you could also use a slope, e.g. when does the area between the linear
    regression and the real data increase by less than say 10%

-effectively, I think we can use a 'loudness' threshold because thats
effectively what a 'listening mother' will actually hear, and therefore is
whats functionally relevant
- we can ignore the contour, it appears that uvseg has a very easy and
intuitive way of grabbing contours.  SO the to do list:

1. reject boxes that dont have enough 'intensity'- this will be a batch
process that modfies the output files that deepsqueak makes
    a. to do this first get a spectrogram of the full day, i think there
    is an option for median, and then use that as a background subtract for
    all pixels.
    b. now get an average peak intensity for all calls in that day, and
    subtract anything less than say 3sd below that, or like less than say
    50% of the distance from that to the mean.
    c. assign all those calls 'reject' calls, or delete those rows

2. regather the contours based on the usvseg method, and then see whether
we need to combine boxes or expand boxes.

3. use uvseg to gather the principal contour(s) and analyze those for
frequency- maybe look at principal frequencie(s) using a findpeaks method,
and then run regressions on those lines to get pitch shifts vs slides

4. use a machine learning algorithm, probably using the full column to
segment usvs, dont fuck with it, it should just work.



%}

%% run this procedure to gather ROIs and remove squeaks that are too soft

% one remaining problem is that I still have low frequency noise being
% detected as signal.  I think i may be able to get over this by running it
% over a gradient smoothing algorithm.  the idea is that there should be a
% high level of at least somewhat locally high gradients that the noise
% wont have.  I can probably get around this by measuring the skew (geomean
% over mean) over say the 10,000 hz surrounding the spline.  the spline
% also shouldnt be too scattered...

edit removeSoftCalls

try
    load('USVmetadata.mat');
catch
    edit USVmetaAggregator 
end

edit mergeMetaCalls

%% just curious about weight
weightmat={}; allweights={};
for i=1:max(USVSession.Cohort)
    uniquedays=unique(USVSession.Age(USVSession.Cohort==i));
    allweights{i}=uniquedays;
    for j=1:length(uniquedays)
        weightmat{i}(:,j)=USVSession.Weight(USVSession.Cohort==i & USVSession.Age==uniquedays(j));
    end
end
figure;
for i=1:length(weightmat)
    errorbar(allweights{i},nanmean(weightmat{i}),nanstd(weightmat{i}));
    hold on;
    keyboard
end
% it appears that cohort 5 is older than reported...
%%

load('C:\Users\John Bladon\Desktop\USVdataFULL2021-10-15.mat');
% this is what deepsqueak gives me to work with...


%% input data for the clustering algorithms
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

For Kmeans:
1. reshapedX
2. zscore slope of 

%}

%% this is the beginning of the Procedure file
%
% produced by zach and jay
%
%% for zachs data heres how they got analyzed:

% zach has a file that he loads first... but he also exported the excel
% data

% this imports the excel calls into a matlab variable
edit importdeepsqueakexcel
% this parses the filenames
edit preprocess_deepsqueak
% this is inside the above to pull data
edit scrape_fileinfo
%%
%%%%%%%%  now analyze data %%%%%%%%%%%%%
edit calls_count_analysis
edit full_plot_analysis