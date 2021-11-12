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

end

