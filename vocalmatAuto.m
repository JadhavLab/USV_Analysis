function [outputArg1,outputArg2] = vocalmatAuto(filename,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here'


% algorithm notes block
%{
Vocalmat uses scripting to run, so it doesnt have a gui and the code
involves a lot of hardcoded parameters and variable deletions/overwrites.

This is an attempt to clean that up and describe/rationalize the algorithm.
It apparently works better than usvseg and deepsqueak on their datasets,
but that isnt to say it works better on mine.  I'll use all three and
decide which is best.

The outline:

1. Check dependencies
    Required libraries: signal processing, image processing, stats and
    machine learning, deep learning

2. Run vocalmat_identifier
    a. load audiofile
    b. parse into one minute sgegments
    c. overlap segments
    d. preallocate F,T,a, grain
    e. for each segment:
        i.   create spectrogram: 1024 nfft, using a 256 pt hamming window, and
             a stepsize of half that(128) pts.  This is not related to time
             though, fo rour data its 2/3 msec, step 1/3 msec. nfft ends up
             being This is a 187 Hz step roughly (linspace(0,nyq,1024)
        ii.  bracket to between 15k and 90k
        iii. convert to psd (10log10)
        iv.  remove half of the padding
        v.   first adjust by subtracting MAX overall, then normalize to 0-1
             Key; This saturates top and bottom 1% of the image and thus
             kills the tails a bit
        vi.  remove last ov the padding
        vii. binarize image using an adaptive sensitivity of .2, and a
             neighborhood of size computed as 2*floor(size(I)/16)+1) so 64
             pix in y direction, and 3.75 seconds in the x direction, low
             high sensitivity of .2 (lets alot of false+ through
       viii. convolve with some shapes: a rectangle of 4x2, then a line of
             4x1: basically its an 'erosion' followed by a 'dilation'...
             first convolve with kernel find local maxima of convolution
             then dilate that by multiplying it by the kernel (to 'open')
         ix. The last few steps are a thresholding, so use the above open
             image as a mask, now you have only your large calls
          x. now gather all calls of 20 pixels or bigger  or...
                (x/.3msec * y/180hz)=20
             average is say 6 pix y, so about 1ms and 1khz bandwidth
         xi. gather the blobs that are at least 4connected (filling holes)
             get their area, threshold AGAIN to 60 pixels (or like 3msec
             and 1khz
    f. Analyze the vocalizations:
        i. gather the x and y coordinates to time freq points
        ii. get the number of frequency jumps
        iii. update local max-min freqs
        iv. get intensity
        v. get each selected freqs (or basicaly harmonic) freq and
        intensity
        vi. local median filtering:

    g. they classify the calls:
        i. this basically means taking the raw image, need to check the
        filtering, and then classifying the raw image based on an
        imagestore.  it is good at finding noise, but thats about it... i
        ran on 4 calls, already misclassified one...

IN SUMMARY:
-they have a comparably sophisticated method to usvseg of detecting calls
- they classify calls when they detect them- into about 12 categories- flat,
up, down, up chev, down chev, hi lo, lo hi, hi lo hi, swoopy, and multiple
hi los, noise and short
-they have a comparably unsophisticated method to deepsqueak of classifying calls
-in classifying calls, we have two options i think.
    1. Find a way to warp the image that makes it LOOK more resolved, then
    train a classifier on some images.
    2. use a set of descriptive statistics to parse the calls, and then
    cluster them in an intelligent way.



%}

% -- save the output from the identifier, in case you only want to rerun the classifier
save_output_files = 0;
% -- max_interval: maximum allowed interval between points to be considered part of one vocalization
max_interval = 20;
% -- minimum_size: minimum number of points to be considered a vocalization
minimum_size = 6;

save_plot_spectrograms    = 0; % plots the spectograms with axis
save_excel_file           = 1; % save output excel file with vocalization stats
scatter_step              = 3; % plot every third point overlapping the vocalization (segmentation)
axes_dots                 = 1; % show the dots overlapping the vocalization (segmentation)
bin_size                  = 300; % in seconds



[a,b,c]=vocalmatIdentifierFx(inputs);


[d,e,f]=vocalmatClassifierFx(inputs2);
end



%
%
%   VOCALMAT IDENTIFIER AND ASSOCIATED FUNCTIONS
%
%
%


function [a,b,c]=vocalmatIdentifierFx(inputs)
raiz = pwd;
disp('[vocalmat]: choose the audio file to be analyzed.');
[vfilename,vpathname] = uigetfile({'*.wav'},'Select the sound track');
cd(vpathname);
p = mfilename('fullpath');

save_spectrogram_background = 0;

% -- local_median: option to use the local median method to detect noise
local_median = 1;


vfile = fullfile(vpathname, vfilename);
mkdir(vfile(1:end-4))
clear time_vocal freq_vocal intens_vocal time_vocal_nogaps freq_vocal_nogaps intens_vocal_nogaps

% -- y1: sampled data; fs: sample rate
[y1,fs] = audioread(vfile);

% -- duration: number of one minute segments in the audio file
duration = ceil(size(y1,1)/(60*fs));
disp(['[vocalmat]: ' vfilename ' has around ' num2str(duration-1) ' minutes.'])

% -- segm_size: duration of each segment to be processed individually, in minutes
% -- overlap: amount of overlap between segments, in seconds
segm_size = 1;
overlap   = 5;
segments  = segm_size:segm_size:duration;
if segments(end) < duration
    segments = [segments, duration];
end

% -- pre-allocate known-size variables for faster performance
num_segments = size(segments,2);
F_orig       = [];
T_orig       = cell(1, num_segments);
A_total      = cell(1, num_segments);
grain_total  = cell(1, num_segments);

for minute_frame = 1:num_segments
% -- run through each segment, compute the spectrogram, and process its outputs

    clear A B y2 S F T P q vocal id grain

    if minute_frame == 1
    % -- y2: current minute frame in seconds, cropped from the whole audio file (y1)
    % -- boundary conditions for first minute, last minute, and files smaller than one minute
        try
            y2 = y1(60*(segments(minute_frame)-segm_size)*fs+1:(60*segments(minute_frame)+overlap)*fs);
        catch
            y2 = y1(60*(segments(minute_frame)-segm_size)*fs+1:end);
        end
    elseif minute_frame == size(segments,2)
        y2 = y1((60*(segments(minute_frame-1))-overlap)*fs+1:end);
    else
        if (60*segments(minute_frame)+overlap)*fs>size(y1,1)
            y2 = y1((60*(segments(minute_frame)-segm_size)-overlap)*fs+1:end);
        else
            y2 = y1((60*(segments(minute_frame)-segm_size)-overlap)*fs+1:(60*segments(minute_frame)+overlap)*fs);
        end
    end

    disp(['[vocalmat][segment (' num2str(minute_frame) ')]: computing the spectrogram.'])
    % -- compute the spectrogram
    % -- nfft: number of points for the Discrete Fourier Transform
    % -- window: windowing function
    % -- nover: number of overlapped samples
    nfft      = 1024; % so 1024 points between fs/2 (180000) and window/2, or 256
    window    = hamming(256); % in my case ~2/3 msec
    nover     = (128); % in my case 1/3 second
    [S,F,T,P] = spectrogram(y2, window, nover, nfft, fs, 'yaxis');
    
    % -- remove frequencies bellow 15kHz
    okFreq  = find(F>15000 & F<90000);

    F = F(okFreq);
    S = S(okFreq,:);
    P = P(okFreq,:);

    % -- convert power spectral density to dB
    P(P==0)=1;
    A = 10*log10(P);
    if minute_frame == 1
        % -- remove first 0.3s of recording (recordings might have abnormal behaviour in this range)
        A = A(:,600:end);
        T = T(:,600:end);
    end
    
    % -- normalize 'A', subtract the maximum value pixel-wise (imcomplement), and adjust contrast (imadjust)
    B = imadjust(imcomplement(abs(A)./max(abs(A(:)))));
    
    % -- adjust minute frame to remove extra padding
    if segments(minute_frame) == segm_size
        F_orig = F;
        lim_inferior = 1;
        lim_superior = find(T<=60*segments(minute_frame),1,'last');
    elseif minute_frame == size(segments,2)
        T = T+(60*(segments(minute_frame-1))-overlap)*ones(size(T,2),1)';
        lim_inferior = find(T>=(60*(segments(minute_frame-1))),1,'first');
        lim_superior = size(T,2);
    else
        T = T+(60*(segments(minute_frame)-segm_size)-overlap)*ones(size(T,2),1)';
        lim_inferior = find(T>=(60*(segments(minute_frame)-segm_size)),1,'first');
        lim_superior = find(T<=60*segments(minute_frame),1,'last');
    end

    T = T(lim_inferior:lim_superior);
    T_orig{minute_frame} = T;
    
    A = A(:,lim_inferior:lim_superior);
    A_total{minute_frame} = A;

    % -- binarize image using an adaptive threshold
    BW = imbinarize(B, 'adaptive', 'Sensitivity', 0.200000, 'ForegroundPolarity', 'bright');
    
    verbose=0;
    if verbose
        figure; subplot(3,1,1); imagesc(T(1:10000),F,A(:,1:10000));
        subplot(3,1,2); imagesc(T(1:10000),F,B(:,1:10000));
        subplot(3,1,3); imagesc(T(1:10000),F,BW(:,1:10000));
    end
    
    % -- morphological image operations
    % -- se: structuring element - rectangle structuring element of size 4x2 pixels
    dimensions = [4 2];
    se = strel('rectangle', dimensions);
    BW = imopen(BW, se);
    
    % -- se: structuring element - line of 4 pixels in length at a 90 degree angle
    length = 4.000000;
    angle  = 90.000000;
    se     = strel('line', length, angle);
    BW     = imdilate(BW, se);
    
    % -- apply mask to original image
    maskedImage      = B;
    maskedImage(~BW) = 0;
    B = maskedImage;
    
    disp(['[vocalmat][segment (' num2str(minute_frame) ')]: computing connected components.'])
    % -- calculate connected components using 4-connected neighborhood policy
    cc = bwconncomp(B, 4);

    % -- calculate area of connected components
    % -- if area is lower than 20, remove
    graindata = regionprops(cc,'Area');
    min_area  = find([graindata.Area]>20) ; % this is hardcoded, need to fix
    grain     = false(size(B));
    for k=1:size(min_area,2)
        grain(cc.PixelIdxList{min_area(k)}) = true;
    end
    grain2 = grain(:,lim_inferior:lim_superior);
    
    disp(['[vocalmat][segment (' num2str(minute_frame) ')]: refining connected components.'])
    % -- recalculate connected components
    % -- if area is lower than 60, remove
    cc        = bwconncomp(grain2, 4);
    graindata = regionprops(cc,'Area');
    
    clear grain grain2;

    min_area  = find([graindata.Area]>60) ;
    grain     = false(size(A));
 
    for k=1:size(min_area,2)
        grain(cc.PixelIdxList{min_area(k)}) = true;
    end
    
    % -- se: line of the 3 pixels in length at a 0 degree angle
    % -- dilate using structuring element 'se'
    length = 3.000000;
    angle  = 0;
    se     = strel('line', length, angle);
    grain  = imdilate(grain, se);
    grain_total{minute_frame} = grain;
end

% -- convert cell array to conventional array
T_orig      = cell2mat(T_orig);
A_total     = cell2mat(A_total);
grain_total = cell2mat(grain_total);

% -- calculate connected components using 4-connected policy, then calculate Area and PixelList of the region
grain       = grain_total;
cc_2        = bwconncomp(grain, 4);
graindata_2 = regionprops(cc_2,'Area','PixelList');


% ----------------------------------------------------------------------------------------------
% -- (3) POST-PROCESSING BEGIN -----------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% -- initialize variables
time_vocal         = [];
id                 = 1;
cc_count           = size(graindata_2,1)-1;
centroid_to_id     = cell(cc_count, 1);

for k = 1:cc_count
% -- for each connected component, get vocalization x-coordinates (time_vocal) and save frequency points (freq_vocal, y-coordinates)
    if k == 1
        time_vocal{id}    = [];
        time_vocal{id}    = unique(graindata_2(k).PixelList(:,1))';
        freq_vocal{id}{1} = [];
        for freq_per_time = 1:size(time_vocal{id},2)
            freq_vocal{id}{freq_per_time} = find(grain(:,time_vocal{id}(freq_per_time))==1);
        end
    else
        if min(graindata_2(k).PixelList(:,1)) - max(time_vocal{id}) > max_interval
        % -- if two points are distant enough, identify as a new vocalization
            id = id + 1;
            time_vocal{id}    = [];
            time_vocal{id}    = unique(graindata_2(k).PixelList(:,1))';
            freq_vocal{id}{1} = [];
            for freq_per_time = 1:size(time_vocal{id},2)
                freq_vocal{id}{freq_per_time} = find(grain(:,time_vocal{id}(freq_per_time))==1);
            end
        else
            time_vocal{id}    = unique([time_vocal{id}, graindata_2(k).PixelList(:,1)']);
            freq_vocal{id}{1} = [];
            for freq_per_time = 1:size(time_vocal{id},2)
                freq_vocal{id}{freq_per_time} = find(grain(:,time_vocal{id}(freq_per_time))==1);
            end
        end
    end
    centroid_to_id{k} = [id, k, T_orig(time_vocal{id}(1)), graindata_2(k).Area];
end

centroid_to_id = cell2mat(centroid_to_id);
centroid_orig  = centroid_to_id;
temp           = [];
idx            = unique(centroid_to_id(:,1));
for k=1:size(idx,1)
    aux  = centroid_to_id((centroid_to_id(:,1)==idx(k)),:);
    temp = [temp; [aux(1,[1 3]) sum(aux(:,4))]];
end
centroid_to_id = temp;

if size(time_vocal,2)>0
% -- if there are vocalizations, remove the ones that have less than 6 points
    disp(['[vocalmat]: removing small vocalizations (less than ' num2str(minimum_size) ' points).'])
    for k=1:size(time_vocal,2)
        if  size(time_vocal{k},2) < minimum_size
            time_vocal{k} = [];
            freq_vocal{k} = [];
        end
    end

    % -- do some cleaning, remove empty cells
    time_vocal = time_vocal(~cellfun('isempty',time_vocal));
    freq_vocal = freq_vocal(~cellfun('isempty',freq_vocal));
    
    freq_harmonic = {};
    time_harmonic = {};

    for k=1:size(time_vocal,2)
    % -- for each vocalization, convert x|y-coordinates to frequency-time points
        for col = 1:size(time_vocal{k},2)
            list_vocal_freq    = find(grain(:,time_vocal{k}(col))==1);
            freq               = F_orig(list_vocal_freq);
            freq_vocal{k}{col} = freq;
            time_vocal{k}(col) = T_orig(time_vocal{k}(col));
        end
    end
        
    for k=1:size(time_vocal,2)
    % -- for each vocalization, check in each timestamp if there is a jump in frequency (harmonic)
        max_local_freq(k) = 0;
        min_local_freq(k) = 200000;
        for time_stamp = 1:size(time_vocal{k},2)
            temp = [];
            if any((freq_vocal{k}{time_stamp} - circshift(freq_vocal{k}{time_stamp} ,[1,0])) > 1000)
            % -- if there are a jumps in frequency, get all jumps
                idx_harmonic = find((freq_vocal{k}{time_stamp} - circshift(freq_vocal{k}{time_stamp} ,[1,0])) > 1000);
                for j=1:size(idx_harmonic,1)
                % -- for each jump, get all frequency points for each range
                    if size(idx_harmonic,1)==1
                    % -- if there's only one jump, get both ranges
                        temp = [temp ; median((freq_vocal{k}{time_stamp}(1:idx_harmonic(j)-1)))];
                        temp = [temp ; median((freq_vocal{k}{time_stamp}(idx_harmonic(j):end)))];
                    else
                    % -- else, sweep jumps and get each range
                        if j==1
                            temp = [temp ; median((freq_vocal{k}{time_stamp}(1:idx_harmonic(j)-1)))];
                        else
                            try
                                temp = [temp ; median((freq_vocal{k}{time_stamp}(idx_harmonic(j-1):idx_harmonic(j)-1)))];
                            catch
                                temp = [temp ; median((freq_vocal{k}{time_stamp}(idx_harmonic(j-1):end)))];
                            end
                        end
                    end
                end

                % -- update local maximum and minimum frequencies
                freq_vocal{k}{time_stamp} = temp;
                if max(temp)>max_local_freq(k)
                    max_local_freq(k) = max(temp);
                end
                if min(temp)<min_local_freq(k)
                    min_local_freq(k) = min(temp);
                end
            else
            % -- if there are no jumps in frequency, only update local maximum and minimum frequencies
                if max((freq_vocal{k}{time_stamp}))>max_local_freq(k)
                    max_local_freq(k) = max((freq_vocal{k}{time_stamp}));
                end
                if min((freq_vocal{k}{time_stamp}))<min_local_freq(k)
                    min_local_freq(k) = min(min((freq_vocal{k}{time_stamp})));
                end
                freq_vocal{k}{time_stamp} = median((freq_vocal{k}{time_stamp}));
            end
        end
    end
    
    for k=1:size(time_vocal,2)
    % -- for each vocalization, get intensity (dB)
        intens_vocal{k} = [];
        for col = 1:size(time_vocal{k},2)
        % -- for each timestamp in a vocalization, select the frequencies belonging to that vocalization
            time_selected                 = time_vocal{k}(col)==T_orig;
            [~, time_selected] = max(time_selected);
            for col2 = 1:size(freq_vocal{k}{col},1)
            % -- for each selected frequency, get its intensity
                freq_selected      = abs(F_orig - freq_vocal{k}{col}(col2));
                [~, freq_selected] = min(freq_selected);
                intens_vocal{k}    = [intens_vocal{k}; A_total(freq_selected,time_selected)];
            end
        end
        
        % -- update local maximum and minimum frequencies
        aux               = abs(F_orig-max_local_freq(k));
        [~, aux]          = min(aux);
        max_local_freq(k) = F_orig(aux);
        
        aux               = abs(F_orig-min_local_freq(k));
        [~, aux]          = min(aux);
        min_local_freq(k) = F_orig(aux);
    end

    median_stats = [];
    
    if local_median == 1
    % -- remove noise using local median
    disp(['[vocalmat]: removing noise by local median.'])
        for k=1:size(time_vocal,2)
        % -- for each vocalization, save timestamp where vocalization begins, connected component area, number of elements, ...
            aux_median_stats = [];
            skip_max_freq    = 0;
            try
               pos              = ceil(size(time_vocal{k},2)/2);
               % pull median db from 200 tstamps before and 200 after call
               % (or roughly 66 msecs?)
               median_db        = median(median(A_total(find(min_local_freq(k)==F_orig)-5:find(max_local_freq(k)==F_orig)+5,...
                   find(T_orig==time_vocal{k}(pos))-200 : find(T_orig==time_vocal{k}(pos))+200)));
               aux_median_stats = [aux_median_stats, time_vocal{k}(1)];
               aux_median_stats = [aux_median_stats, centroid_to_id(find(centroid_to_id(:,2)==time_vocal{k}(1)),3)];
               % spal frequency up and down 5 bins (or about 1000 hz)
               aux_median_stats = [aux_median_stats, numel(A_total(find(min_local_freq(k)==F_orig)-5:find(max_local_freq(k)==F_orig)+5,find(T_orig==time_vocal{k}(pos))-200 : find(T_orig==time_vocal{k}(pos))+200))];
            catch
            % -- boundary conditions
                pos = ceil(size(time_vocal{k},2)/2);
                if find(min_local_freq(k)==F_orig)-5 < 1
                % -- check frequency point is not out or range (lower bound)
                    if find(max_local_freq(k)==F_orig)+5 > size(A_total,1)
                        okFreq = 1;
                        max_freq = size(F_orig,1);
                        max_time = find(T_orig==time_vocal{k}(pos))+200;
                        min_time = find(T_orig==time_vocal{k}(pos))-200;
                        if min_time<1
                            min_time=1;
                        end
                    else
                        okFreq = 1;
                        max_freq = find(max_local_freq(k)==F_orig)+5;
                        max_time = find(T_orig==time_vocal{k}(pos))+200;
                        min_time = find(T_orig==time_vocal{k}(pos))-200;
                        if min_time<1
                            min_time=1;
                        end
                    end
                    skip_max_freq = 1;
                end
                if find(max_local_freq(k)==F_orig)+5 >= size(A_total,1) && skip_max_freq==0
                % -- check frequency point is not out or range (upper bound)
                    max_freq = size(F_orig,1);
                    okFreq = find(min_local_freq(k)==F_orig)-5;
                    max_time = find(T_orig==time_vocal{k}(pos))+200;
                    min_time = find(T_orig==time_vocal{k}(pos))-200;
                    if min_time < 1
                        min_time=1;
                    end
                end
                if find(T_orig==time_vocal{k}(pos))-200 < 1
                % -- check time point is not out or range (lower bound)
                    min_time = 1;
                    max_time = find(T_orig==time_vocal{k}(pos))+200;
                    max_freq = find(max_local_freq(k)==F_orig)+5;
                    if max_freq > size(A_total,1)
                        max_freq = size(A_total,1);
                    end
                    okFreq = find(min_local_freq(k)==F_orig)-5;
                    if okFreq < 1
                        okFreq = 1;
                    end
                end
                if find(T_orig==time_vocal{k}(pos))+200 >= size(A_total,2)
                % -- check time point is not out or range (upper bound)
                    max_time = size(A_total,2);
                    min_time = find(T_orig==time_vocal{k}(pos))-200;
                    max_freq = find(max_local_freq(k)==F_orig)+5;
                    if max_freq > size(A_total,1)
                        max_freq = size(A_total,1);
                    end
                    okFreq = find(min_local_freq(k)==F_orig)-5;
                    if okFreq < 1
                        okFreq=1;
                    end
                end
                median_db        = median(median(A_total(okFreq:max_freq,min_time:max_time)));
                aux_median_stats = [aux_median_stats, time_vocal{k}(1)];
                aux_median_stats = [aux_median_stats, centroid_to_id(find(centroid_to_id(:,2)==time_vocal{k}(1)),3)];
                aux_median_stats = [aux_median_stats, numel(A_total(okFreq:max_freq,min_time:max_time))];
                
            end

            temp             = sort(intens_vocal{k});
            aux_median_stats = [aux_median_stats, size(temp,1)];
            aux_median_stats = [aux_median_stats, [median(temp(end-5:end))  median_db]];
            elim_by_median   = 0;
            aux_median_stats = [aux_median_stats, elim_by_median];
            median_stats(k,:) = aux_median_stats;
        end
        % median stats is... 1
        ratio = median_stats(:,5)./median_stats(:,6);
        [y,t]=ecdf(ratio);
        aux = round(linspace(1,size(t,1),35)); % Downsample to 50 points only
        t = t(aux);
        y = y(aux);
        K=LineCurvature2D([t,y]);
        K = K*10^-3;
        [maxx maxx] = max(K);
        th_ratio = t(maxx);

        if th_ratio<0.9
            th_ratio=0.92;
        end
        disp(['[vocalmat]: minimal ratio = ' num2str(th_ratio) '.'])
        
        for k=1:size(time_vocal,2)
            if median_stats(k,5) < th_ratio*median_stats(k,6)
                time_vocal{k}=[];
                freq_vocal{k}=[];
                intens_vocal{k}=[];
                median_stats(k,7) = 1;
            end
        end
        
        % -- do some cleaning, remove empty cells
        time_vocal        = time_vocal(~cellfun('isempty',time_vocal));
        freq_vocal        = freq_vocal(~cellfun('isempty',freq_vocal));
        intens_vocal      = intens_vocal(~cellfun('isempty',intens_vocal));
        intens_vocal_orig = intens_vocal;
        freq_vocal_orig   = freq_vocal;
        time_vocal_orig   = time_vocal;
        
    end
    intens_orig = intens_vocal;
    
    for k=1:size(time_vocal,2)
        temp={};
        for kk=1:size(freq_vocal{k},2)
        % -- for each vocalization, order its intensities in the same pattern as its frequencies
            temp = [ temp intens_vocal{k}(1:size(freq_vocal{1,k}{1,kk},1))];
            intens_vocal{k}(1:size(freq_vocal{1,k}{1,kk},1)) = [];
        end
        intens_vocal{k} = temp;
    end

    vfilename  = vfilename(1:end-4);
    if save_output_files == 1
        disp(['[vocalmat]: saving output files.'])
        % -- output identified vocalizations
%         cd(fullfile(root_path, 'audios'))
        
        save(fullfile(vfile(1:end-4), ['output_short_' vfilename]), 'T_orig', 'F_orig', 'time_vocal', 'freq_vocal', 'vfilename', 'intens_vocal', 'median_stats')
        
        if save_spectrogram_background == 1
            save(fullfile(vfile(1:end-4), ['output_' vfilename]), 'T_orig', 'F_orig', 'time_vocal', 'freq_vocal', 'vfilename', 'intens_vocal', 'median_stats', 'A_total', '-v7.3', '-nocompression')
        end
    end
    
    warning('off', 'MATLAB:save:sizeTooBigForMATFile')
    clear y y1 S F T P fs q nd vocal id
    
    toc
end
disp(['[vocalmat]: ' vfilename ' has ' num2str(size(time_vocal,2)) ' vocalizations.'])
end


function k=LineCurvature2D(Vertices,Lines)
% This function calculates the curvature of a 2D line. It first fits 
% polygons to the points. Then calculates the analytical curvature from
% the polygons;

%  k = LineCurvature2D(Vertices,Lines)
% 
% inputs,
%   Vertices : A M x 2 list of line points.
%   (optional)
%   Lines : A N x 2 list of line pieces, by indices of the vertices
%         (if not set assume Lines=[1 2; 2 3 ; ... ; M-1 M])
%
% outputs,
%   k : M x 1 Curvature values
%
% Example, Circle
%  r=sort(rand(15,1))*2*pi;
%  Vertices=[sin(r) cos(r)]*10;
%  Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
%  k=LineCurvature2D(Vertices,Lines);
%
%  figure,  hold on;
%  N=LineNormals2D(Vertices,Lines);
%  k=k*100;
%  plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
%  plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
%  plot(sin(0:0.01:2*pi)*10,cos(0:0.01:2*pi)*10,'r.');
%  axis equal;
%
% Example, Hand
%  load('testdata');
%  k=LineCurvature2D(Vertices,Lines);
%
%  figure,  hold on;
%  N=LineNormals2D(Vertices,Lines);
%  k=k*100;
%  plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
%  plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
%  plot(Vertices(:,1),Vertices(:,2),'r.');
%  axis equal;
%
% Function is written by D.Kroon University of Twente (August 2011)

% If no line-indices, assume a x(1) connected with x(2), x(3) with x(4) ...
if(nargin<2)
    Lines=[(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
end

% Get left and right neighbor of each points
Na=zeros(size(Vertices,1),1); Nb=zeros(size(Vertices,1),1);
Na(Lines(:,1))=Lines(:,2); Nb(Lines(:,2))=Lines(:,1);

% Check for end of line points, without a left or right neighbor
checkNa=Na==0; checkNb=Nb==0;
Naa=Na; Nbb=Nb;
Naa(checkNa)=find(checkNa); Nbb(checkNb)=find(checkNb);

% If no left neighbor use two right neighbors, and the same for right... 
Na(checkNa)=Nbb(Nbb(checkNa)); Nb(checkNb)=Naa(Naa(checkNb));

% Correct for sampeling differences
Ta=-sqrt(sum((Vertices-Vertices(Na,:)).^2,2));
Tb=sqrt(sum((Vertices-Vertices(Nb,:)).^2,2)); 

% If no left neighbor use two right neighbors, and the same for right... 
Ta(checkNa)=-Ta(checkNa); Tb(checkNb)=-Tb(checkNb);

% Fit a polygons to the vertices 
% x=a(3)*t^2 + a(2)*t + a(1) 
% y=b(3)*t^2 + b(2)*t + b(1) 
% we know the x,y of every vertice and set t=0 for the vertices, and
% t=Ta for left vertices, and t=Tb for right vertices,  
x = [Vertices(Na,1) Vertices(:,1) Vertices(Nb,1)];
y = [Vertices(Na,2) Vertices(:,2) Vertices(Nb,2)];
M = [ones(size(Tb)) -Ta Ta.^2 ones(size(Tb)) zeros(size(Tb)) zeros(size(Tb)) ones(size(Tb)) -Tb Tb.^2];
invM=inverse3(M);
a(:,1)=invM(:,1,1).*x(:,1)+invM(:,2,1).*x(:,2)+invM(:,3,1).*x(:,3);
a(:,2)=invM(:,1,2).*x(:,1)+invM(:,2,2).*x(:,2)+invM(:,3,2).*x(:,3);
a(:,3)=invM(:,1,3).*x(:,1)+invM(:,2,3).*x(:,2)+invM(:,3,3).*x(:,3);
b(:,1)=invM(:,1,1).*y(:,1)+invM(:,2,1).*y(:,2)+invM(:,3,1).*y(:,3);
b(:,2)=invM(:,1,2).*y(:,1)+invM(:,2,2).*y(:,2)+invM(:,3,2).*y(:,3);
b(:,3)=invM(:,1,3).*y(:,1)+invM(:,2,3).*y(:,2)+invM(:,3,3).*y(:,3);

% Calculate the curvature from the fitted polygon
k = 2*(a(:,2).*b(:,3)-a(:,3).*b(:,2)) ./ ((a(:,2).^2+b(:,2).^2).^(3/2));
end



function  Minv = inverse3(M)
% This function does inv(M) , but then for an array of 3x3 matrices
adjM(:,1,1)=  M(:,5).*M(:,9)-M(:,8).*M(:,6);
adjM(:,1,2)=  -(M(:,4).*M(:,9)-M(:,7).*M(:,6));
adjM(:,1,3)=  M(:,4).*M(:,8)-M(:,7).*M(:,5);
adjM(:,2,1)=  -(M(:,2).*M(:,9)-M(:,8).*M(:,3));
adjM(:,2,2)=  M(:,1).*M(:,9)-M(:,7).*M(:,3);
adjM(:,2,3)=  -(M(:,1).*M(:,8)-M(:,7).*M(:,2));
adjM(:,3,1)=  M(:,2).*M(:,6)-M(:,5).*M(:,3);
adjM(:,3,2)=  -(M(:,1).*M(:,6)-M(:,4).*M(:,3));
adjM(:,3,3)=  M(:,1).*M(:,5)-M(:,4).*M(:,2);
detM=M(:,1).*M(:,5).*M(:,9)-M(:,1).*M(:,8).*M(:,6)-M(:,4).*M(:,2).*M(:,9)+M(:,4).*M(:,8).*M(:,3)+M(:,7).*M(:,2).*M(:,6)-M(:,7).*M(:,5).*M(:,3);
Minv=bsxfun(@rdivide,adjM,detM);
end




%
%
%   VOCALMAT CLASSIFIER AND ASSOCIATED FUNCTIONS
%
%
%

function [d,e,f]=vocalmatClassifierFx(inputs2)


size_spectrogram   = [227 227];
use_DL             = 1;
plot_stats_per_bin = 1;


raiz = pwd;
% find the preconfigured classifier model
model_class_DL = load('Mdl_categorical_DL.mat');
model_class_DL = model_class_DL.netTransfer;

% [vfilename,vpathname] = uigetfile({'*.mat'},'Select the output file')
% disp(['Reading ' vfilename])
vfile = fullfile(vpathname,vfilename); 
% load(vfile);
%cd(vpathname);
%list = dir('*output*.mat');
%diary(['Summary_classifier' num2str(horzcat(fix(clock))) '.txt'])

%Setting up
p = mfilename('fullpath');
fprintf('\n')


%We are gonna get only 10 points (time stamps) to classify the vocalization

% Grimsley, Jasmine, Marie Gadziola, and Jeff James Wenstrup. 
% "Automated classification of mouse pup isolation syllables: 
% from cluster analysis to an Excel-based mouse pup syllable 
% classification calculator." Frontiers in behavioral neuroscience 6 (2013): 89.

%     disp('Verify vocalizations for steps')
stepup_count=[];
stepdown_count=[];
harmonic_count=[];
flat_count=[];
chevron_count=[];
revchevron_count=[];
downfm_count=[];
upfm_count=[];
complex_count=[];
noisy_vocal_count=[];
nonlinear_count = [];
short_count = [];
noise_count = [];
noise_count_dist = [];
corr_yy2_yy3 = [];
corr_yy2_yy4 = [];
max_prom = [];
max_prom2 = [];
duration = [];


% disp('[vocalmat][classifier]: checking for empty cells')
time_vocal = time_vocal(~cellfun('isempty',time_vocal));
freq_vocal = freq_vocal(~cellfun('isempty',freq_vocal));
intens_vocal = intens_vocal(~cellfun('isempty',intens_vocal));

output=[];
cd(vpathname)
if ~exist(vfilename, 'dir')
    mkdir(vfilename)
end
cd(vfilename)

disp('[vocalmat][classifier]: running analysis!')

for k=1:size(time_vocal,2)
    
    harmonics = cell(1,size(time_vocal,2));
    
    current_freq = [];
    harmonic_candidate = [];
    skip_current = 0;
    for time_stamp = 1:size(time_vocal{k},2)-1
        
        if size(freq_vocal{k}{time_stamp+1},1)>1 %Probably we have an harmonic
            if (size(freq_vocal{k}{time_stamp},1)>1) %Check if they have same size (could be the continuation of harmonic)
                if time_stamp==1 %If the vocalization starts with an harmonic
                    current_freq = freq_vocal{k}{time_stamp}(1);
                    harmonic_candidate = freq_vocal{k}{time_stamp}(2);
                    if size(harmonic_candidate,1)==1
                        start_harmonic = time_vocal{k}(time_stamp);
                    end
                else
                    aux = freq_vocal{k}{time_stamp+1} - current_freq(end)*ones(size(freq_vocal{k}{time_stamp+1},1),1);
                    [mini,mini]=min(abs(aux));
                    temp = freq_vocal{k}{time_stamp+1};
                    current_freq = [current_freq; temp(mini)]; temp(mini) = [];
                    if size(harmonic_candidate,1)>1
                        if abs(temp - harmonic_candidate(end)) < 10000
                            harmonic_candidate = [harmonic_candidate; temp(1)];
                        else %if it is >10khz then it is already another harmonic
                            if size(harmonic_candidate,1)>10
                                harmonic_count = [harmonic_count;k];
                            end
                            harmonic_candidate = temp;
                        end
                    else
                        harmonic_candidate = [harmonic_candidate; temp(1)];
                    end
                    if size(harmonic_candidate,1)==1
                        start_harmonic = time_vocal{k}(time_stamp);
                    end
                end
            else %Find the closests freq to be the current and classify the other as harmonic candidate
                try
                    aux = freq_vocal{k}{time_stamp+1} - current_freq(end)*ones(size(freq_vocal{k}{time_stamp+1},1),1);
                catch
                    aux = freq_vocal{k}{time_stamp+1} - freq_vocal{k}{time_stamp}*ones(size(freq_vocal{k}{time_stamp+1},1),1);
                end
                
                [mini,mini]=min(abs(aux));
                temp = freq_vocal{k}{time_stamp+1};
                current_freq = [current_freq; temp(mini)]; temp(mini) = [];
                harmonic_candidate = [harmonic_candidate; temp];
                if size(harmonic_candidate,1)==1 || (size(harmonic_candidate,1)>1 && time_stamp==1)
                    start_harmonic = time_vocal{k}(time_stamp);
                end
            end
            
        else %There is nothing similar to harmonic right now... but there was before?
            if (size(freq_vocal{k}{time_stamp},1)>1)
                %                So... Was it an harmonic or not?
                if time_stamp == 1 %If the vocalization starts with something that reminds a vocalization
                    aux = freq_vocal{k}{time_stamp} - freq_vocal{k}{time_stamp+1}*ones(size(freq_vocal{k}{time_stamp},1),1);
                    [mini,mini]=min(abs(aux));
                    temp = freq_vocal{k}{time_stamp};
                    current_freq = [current_freq; temp(mini)]; temp(mini) = [];
                    harmonic_candidate = [harmonic_candidate; temp];
                    if size(harmonic_candidate,1)==1
                        start_harmonic = time_vocal{k}(time_stamp);
                    end
                end
                
                if abs(freq_vocal{k}{time_stamp+1} - harmonic_candidate(end)) < abs(freq_vocal{k}{time_stamp+1} - current_freq(end)) %Continued on the line that we thought was harmonic. So it is not harmonic
                    if size(harmonic_candidate,1)> size(current_freq,1)
                        
                        current_freq = [current_freq; freq_vocal{k}{time_stamp+1}];
                        harmonic_candidate = [];
                    else %current_freq > harmonic_candidate -> So it is a jump, not a harmonic
                        if size(harmonic_candidate,1)>10% && size(harmonic_candidate,1)/ size(current_freq,1)>0.8 %If the harmonic is big and close to the size of current_freq
                            
                            if (time_stamp+2 < size(time_vocal{k},2)) && any(abs(freq_vocal{k}{time_stamp+2} - current_freq(end)) < abs(freq_vocal{k}{time_stamp+2} - harmonic_candidate(end))) %Is there any chance of continuing with the current_freq?
                                harmonic_candidate = [harmonic_candidate; freq_vocal{k}{time_stamp+1}];
                                skip_current = 1;
                                harmonic_count = [harmonic_count;k];
                            else
                                current_freq(end-size(harmonic_candidate,1)+1:end) = harmonic_candidate;
                                current_freq = [current_freq; freq_vocal{k}{time_stamp+1}];
                                harmonic_candidate = [];
                                harmonic_count = [harmonic_count;k];
                            end
                            
                        else %So they just overlapped for a little while, but was actually a step
                            harmonic_candidate = [];
                        end
                    end
                else %It was an harmonic after all
                    current_freq = [current_freq; freq_vocal{k}{time_stamp+1}];
                    if size(harmonic_candidate,1)>10 % at least 10 points to say it was really an harmonic
                        harmonic_count = [harmonic_count;k];
                    end
                    harmonic_candidate = [];
                end
                
            else
                aux = freq_vocal{k}{time_stamp+1} - freq_vocal{k}{time_stamp};
                if skip_current==0
                    current_freq = [current_freq; freq_vocal{k}{time_stamp}];
                end
                skip_current = 0;
                
            end
        end
        
    end
    
    %Extra filtering by removing the points with intensity below 5% of the average
    tabela = [];
    for kk = 1:size(time_vocal{k},2)
        for ll = 1:size(freq_vocal{k}{kk},1)
            tabela = [tabela; time_vocal{k}(kk) freq_vocal{k}{kk}(ll) intens_vocal{k}{kk}(ll)];
        end
    end
    tabela_all_points{k} = tabela;
end

cd(raiz)

if use_DL==1
    if save_plot_spectrograms==1
        fig = figure('Name',vfilename,'NumberTitle','off','Position',[300 200 1167 875]);
    end

    cd(vpathname)
    if ~exist(vfilename, 'dir')
        mkdir(vfilename)
    end
    cd(vfilename)
    
    if (~exist([vfile '\All_axes'],'dir') && save_plot_spectrograms==1)
        mkdir('All_axes')
    end
    
    if ~exist([vfile '\All'],'dir')
        mkdir('All')
    end

    for id_vocal = 1:size(time_vocal,2)
        %         cd(raiz)
        dx = 0.22;
        
        T_min_max = [-dx/2 dx/2]+[time_vocal{id_vocal}(ceil(size(time_vocal{id_vocal},2)/2)) time_vocal{id_vocal}(ceil(size(time_vocal{id_vocal},2)/2))];
        [T_min T_min] = min(abs(T_orig - T_min_max(1)));
        [T_max T_max] = min(abs(T_orig - T_min_max(2)));
        
        if save_plot_spectrograms==1
           if save_plot_spectrograms==1
            clf('reset');
            hold on;
            surf(T_orig(T_min:T_max),F_orig,A_total(:,T_min:T_max),'edgecolor','none');
            axis tight; view(0,90);
            colormap(gray);
            xlabel('Time (s)'); ylabel('Freq (Hz)');
            
            if axes_dots == 1
                for time_stamp = 1:scatter_step:size(time_vocal{id_vocal},2)
                    try
                        scatter(time_vocal{id_vocal}(time_stamp)*ones(size(freq_vocal{id_vocal}{time_stamp}')),freq_vocal{id_vocal}{time_stamp}',[],'b');
                    catch
                        scatter(time_vocal{id_vocal}(time_stamp-1)*ones(size(freq_vocal{id_vocal}{time_stamp-1}')),freq_vocal{id_vocal}{time_stamp}',[],'b');
                    end
                end
            end
            set(gca,'fontsize', 18);
            frame = getframe(fig);
            imwrite(frame.cdata, fullfile(vpathname , vfilename, 'All_axes', [num2str(id_vocal)  '.png']), 'png');
            hold off;
            
            end
        end        
        img = imresize(flipud(mat2gray(A_total(:,T_min:T_max))),size_spectrogram);
        img = cat(3, img, img, img); % why are we making a triplet of this image...
        % need to register this first...
        % img = cat(3,rescale(img(:,:,1)), flipud(filter(-cos(linspace(0,pi,15)),1,flipud(img(:,:,1)))),...
        %    filter(-cos(linspace(0,pi,15)),1,(img(:,:,1))));
        % img(1:15,:,3)=0; img(size(img,1)-15:size(img,1),:,2)=0;
        % i think you would have more power if you filtered the freq domain
        % with an expected bandwith of your usv. in this case looks like
        % cos(linspace(0,pi,15))
        
        % they also dont binarize the image before doing recognition...wtf
        
        imwrite(img,fullfile(vpathname, vfilename, 'All', [num2str(id_vocal)  '.png']))
        
    end
    
end

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rebuild curr_freq
dist_between_points=[];
slopes =[];
all_jumps =[];
higher_jumps=[];
lower_jumps=[];
table_out={};
duration=[];
num_points = 15;
for k=1:size(freq_vocal,2)
    curr_freq = [];
    intens_freq =[];
    for kk=1:size(freq_vocal{k},2)
        if size(freq_vocal{k}{kk},1)>1
            if isempty(curr_freq)
                j=2;
                try
                    while size(freq_vocal{k}{j},1)>1 %find the first element without "harmonic"
                        j=j+1;
                    end
                    [min_idx, min_idx] = min(abs((freq_vocal{k}{j}*ones(size(freq_vocal{k}{kk},1),1)) - freq_vocal{k}{kk}));
                    curr_freq = [curr_freq; freq_vocal{k}{kk}(min_idx)];
                    intens_freq = [intens_freq; intens_vocal{k}{kk}(min_idx)];
                catch
                    curr_freq = [curr_freq; freq_vocal{k}{kk}(1)]; %Just get any element if there is no time stamp with only one
                    intens_freq = [intens_freq; intens_vocal{k}{kk}(1)];
                end
                
            else
                [min_idx, min_idx] = min(abs((curr_freq(end)*ones(size(freq_vocal{k}{kk},1),1)) - freq_vocal{k}{kk}));
                curr_freq = [curr_freq; freq_vocal{k}{kk}(min_idx)];
                intens_freq = [intens_freq; intens_vocal{k}{kk}(min_idx)];
            end
        else
            curr_freq = [curr_freq; freq_vocal{k}{kk}];
            intens_freq = [intens_freq; intens_vocal{k}{kk}];
        end
    end
    curr_freq_total{k} = curr_freq;
end

if use_DL==1
    
    % in 4 images it misclassified one, it classified a complex call as a
    % flat call... this may be only good for detecting noise...
    validationImages = imageDatastore(fullfile(vpathname, vfilename, 'All'));
    [predictedLabels, scores] = classify(model_class_DL,validationImages);
    lista = [validationImages.Files, predictedLabels];
    
    AA2 = cellstr(lista);
    AA = array2table(AA2);
    ttt = model_class_DL.Layers(25).ClassNames;
    ttt2 = cellstr(num2str(2*ones(12,1)));
    s = strcat(ttt,ttt2);
    T2 = array2table(scores,'VariableNames',s');
    
    % AA2 = strsplit(cell2mat(AA2(1,1)),'\');
    for k=1:size(AA2,1)
        AA1 = strsplit(cell2mat(AA2(k,1)),{'/','\'});
        AA3 = str2double(AA1{end}(1:end-4));
        %     AA4 = str2double(AA1{end}(1:end-20));
        AA2(k,3) = num2cell(AA3);
    end
    
    T_classProb = [T2, AA, array2table(cell2mat(AA2(:,3)))];
    T_classProb.Properties.VariableNames{15} = 'NumVocal';
    T_classProb.Properties.VariableNames{14} = 'DL_out';
    T_classProb = sortrows(T_classProb,'NumVocal','ascend');
end

if use_DL==1
%     temp = [T_classProb];
    writetable(T_classProb,fullfile(vfile,[vfilename '_DL.xlsx']))
end
save T_classProb T_classProb
% 
chevron_count = sum(strcmp(T_classProb.DL_out,'chevron'));
complex_count = sum(strcmp(T_classProb.DL_out,'complex'));
down_fm_count = sum(strcmp(T_classProb.DL_out,'down_fm'));
flat_count = sum(strcmp(T_classProb.DL_out,'flat'));
mult_steps_count = sum(strcmp(T_classProb.DL_out,'mult_steps'));
noise_count = sum(strcmp(T_classProb.DL_out,'noise_dist'));
rev_chevron_count = sum(strcmp(T_classProb.DL_out,'rev_chevron'));
short_count = sum(strcmp(T_classProb.DL_out,'short'));
step_down_count = sum(strcmp(T_classProb.DL_out,'step_down'));
step_up_count = sum(strcmp(T_classProb.DL_out,'step_up'));
two_steps_count = sum(strcmp(T_classProb.DL_out,'two_steps'));
up_fm_count = sum(strcmp(T_classProb.DL_out,'up_fm'));
noise_dist_count = sum(strcmp(T_classProb.DL_out,'noise_dist'));
harmonic_count = unique(harmonic_count);
noisy_vocal_count = unique(noisy_vocal_count);

disp(['[vocalmat][classifier]: total number of vocalizations: ' num2str(size(time_vocal,2)-noise_dist_count) ' vocalizations (' num2str(noise_dist_count) ' were noise)']);

for j=1:size(model_class_DL.Layers(25,1).ClassNames)
    eval(['disp([''' cell2mat(model_class_DL.Layers(25,1).ClassNames(j)) ': '' num2str('  cell2mat(model_class_DL.Layers(25,1).ClassNames(j)) '_count)])'])
end

% Fixed up to here.
if save_excel_file==1
    %     names2 = model_class_DL_RF.ClassNames;
    names = [{'Names_vocal'};{'Start_time'}; {'End_time'}; {'Inter_vocal_interval'}; {'Inter_real_vocal_interval'}; {'Duration'}; {'min_freq_main'}; {'max_freq_main'};{'mean_freq_main'};{'Bandwidth'};{'min_freq_total'};...
        {'max_freq_total'};{'mean_freq_total'};{'min_intens_total'};{'max_intens_total'}; {'corrected_max_intens_total'};{'Background_intens'};{'mean_intens_total'};{'Class'};{'Harmonic'};{'Noisy'}];
    tabela = zeros(size(T_classProb,1),size(names,1));
    tabela(:,1) = 1:size(T_classProb,1);
    tabela = num2cell(tabela);
    
    if ~isempty(noisy_vocal_count)
        tabela(noisy_vocal_count,21)= {1};
    end
    
    if ~isempty(harmonic_count)
        tabela(harmonic_count,20)= {1};
    end
    
    for i=1:size(time_vocal,2)
        time_start(i) = time_vocal{i}(1);
        time_end(i) = time_vocal{i}(end);
        if i>1
            time_interval(i) = time_start(i)-time_end(i-1);
        else
            time_interval(i) = NaN;
        end
        duration(i) = time_end(i)-time_start(i);
        if ~isempty(curr_freq_total{i}), min_freq_main(i) = min(curr_freq_total{i}); else min_freq_main(i) = NaN; end
        if ~isempty(curr_freq_total{i}), max_freq_main(i) = max(curr_freq_total{i}); else max_freq_main(i) = NaN; end
        mean_freq_main(i) = mean(curr_freq_total{i});
        min_freq_total(i) = min(tabela_all_points{i}(:,2));
        max_freq_total(i) = max(tabela_all_points{i}(:,2));
        mean_freq_total(i) = mean(tabela_all_points{i}(:,2));
        min_intens_total(i) = min(tabela_all_points{i}(:,3));
        max_intens_total(i) = max(tabela_all_points{i}(:,3));
        mean_intens_total(i) = mean(tabela_all_points{i}(:,3));
    end
    
    tabela(:,19) = T_classProb.DL_out;
    
    noise_idx = strcmp(tabela(:,18),'noise_dist');
    time_start_real = time_start; time_start_real(noise_idx) = NaN;
    time_end_real = time_end; time_end_real(noise_idx) = NaN;
    curr_time = NaN;
    for i=1:size(time_start_real,2)
        if ~isnan(time_start_real(i))
            time_interval_real(i) = time_start_real(i) - curr_time;
            curr_time = time_end_real(i);
        else
            time_interval_real(i) = NaN;
        end
    end
    
    median_stats = [ zeros(size(median_stats,1),1) median_stats, zeros(size(median_stats,1),1)];
    for k=1:size(time_start,2)
        median_stats(find(median_stats(:,2)==time_start(k)),end) = 1;
        median_stats(find(median_stats(:,2)==time_start(k)),1) = k;
    end
    
    median_stats(:,7) = median_stats(:,7)/0.9;
    median_stats = median_stats(median_stats(:,1)>0,:);
    
    
    tabela(:,2) = num2cell(time_start');
    tabela(:,3) = num2cell(time_end');
    tabela(:,4) = num2cell(time_interval');
    tabela(:,5) = num2cell(time_interval_real');
    tabela(:,6) = num2cell(duration');
    tabela(:,7) = num2cell(min_freq_main');
    tabela(:,8) = num2cell(max_freq_main');
    tabela(:,9) = num2cell(mean_freq_main');
    tabela(:,10) = num2cell(max_freq_main'-min_freq_main');
    tabela(:,11) = num2cell(min_freq_total');
    tabela(:,12) = num2cell(max_freq_total');
    tabela(:,13) = num2cell(mean_freq_total');
    tabela(:,14) = num2cell(min_intens_total');
    tabela(:,15) = num2cell(max_intens_total');
    corrected_max_intens_total = max_intens_total' - median_stats(:,7);
    tabela(:,16) = num2cell(corrected_max_intens_total);
    tabela(:,17) = num2cell(median_stats(:,7)');
    tabela(:,18) = num2cell(mean_intens_total');
    
    names = transpose(names);
    T = array2table(tabela);
    T.Properties.VariableNames = names;
    %     VM1_out.Properties.VariableNames{1} = 'VM1_out';
    if exist([vfilename '.xlsx'])>0
        delete([vfilename '.xlsx'])
    end
    
    writetable(T,fullfile(vfile, [vfilename '.xlsx']))
end

% Estimate number of bins given the bin size
aux = ~strcmp(T.Class,'noise_dist');
T_no_noise = T(aux,:);
if size(T_no_noise,1)>0
    num_of_bins = ceil(max(cell2mat(T_no_noise.Start_time))/bin_size);
    edges = 0:bin_size:num_of_bins*bin_size;
    [num_vocals_in_bin] = histcounts(cell2mat(T_no_noise.Start_time),edges);
    
    
    
    disp(['[vocalmat][classifier]: vocalizations per bin (not considering noise):'])
    for k=1:num_of_bins
        disp(['Bin_' num2str(k) '(' num2str(edges(k)) '-' num2str(edges(k+1)) 's): ' num2str(num_vocals_in_bin(k))])
    end
    
    if plot_stats_per_bin ==1
        
        %Show classes per bin
        for j=1:size(model_class_DL.Layers(25, 1).ClassNames  )
            aux = strcmp(T.Class,model_class_DL.Layers(25, 1).ClassNames (j));
            T_class = T(aux,:);
            [num_vocals_in_bin,~] = histcounts(cell2mat(T_class.Start_time),edges);
            disp(['[vocalmat][classifier]: vocalizations per bin for class ' cell2mat(model_class_DL.Layers(25, 1).ClassNames(j)) ' :'])
            for k=1:num_of_bins
                disp(['Bin_' num2str(k) '(' num2str(edges(k)) '-' num2str(edges(k+1)) 's): ' num2str(num_vocals_in_bin(k))])
            end
        end
        
    end
    
else
    disp('[vocalmat][classifier]: no real vocalizations detected in this file.')
end


end