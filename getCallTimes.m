% getCallTimes
% Get Call times for a list of wav files:
% 

% wrapper for detect function

% main function:

% edit usvsegDetect

% and to implement the usvsegDetect:

wavfiledir='G:\USV data\Raw\wav';
wavfilenames=getAllFiles(wavfiledir);
fileinfo=dir(wavfiledir); fileinfo([fileinfo.isdir])=[];
outdir='G:\USV data\segData';
if ~isfolder(outdir)
    mkdir('G:\USV data\segData');
end
nameadd='_callData';
bigclock=tic;
wb=waitbar(0,'Starting to build each session');
for i=427:length(wavfilenames)
    smallclock=tic;
    audiodata=audioinfo(wavfilenames{i});
    
    sessdata=rmfield(fileinfo(i),{'bytes','isdir'});
    myfilename=sessdata.name;
    unders=find(myfilename=='_' | myfilename=='-');
    sessdata.Cohort=str2double(myfilename(2:unders(1)-1)); % cohort
    sessdata.Age=str2double(myfilename(unders(1)+2:unders(2)-1)); % age
    % pull rat name from file name
    sessdata.Ratname=myfilename(unders(2)+1:end-4);
    
    [spect,thrshd,params,onoffset,onoffsetm,blobs]=usvsegDetect(audiodata);
    % spect is the flattened centered im, thrshd is the pix that pass
    % threshold after masking
    
    segCalls=table(onoffset(:,1),onoffset(:,2),ones(size(onoffset,1),1),...
        'VariableNames',{'onsetTime','offsetTime','Accept'});
    
    calldata=getCallStats(segCalls,params,spect,thrshd);
    
    % reduce the data size of the spect
    spect=uint8(rescale(spect)*256);
    save(fullfile(outdir,[sessdata.name(1:end-4) nameadd]),...
        'audiodata','params','spect','blobs','segCalls','calldata');
    waitbar(i/length(wavfilenames),wb,sprintf('Its taken %d minutes, Likely %d more',...
        round(toc(bigclock)/60),round(toc(smallclock)*(length(wavfilenames)-i)/60)));
end
