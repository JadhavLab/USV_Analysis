% getCallTimes
% Get Call times for a list of wav files:
% 

% wrapper for detect function

% main function:

% edit usvsegDetect

% and to implement the usvsegDetect:

% for work machinie
%wavfiledir='G:\USV data\Raw\wav';
%outdir='G:\USV data\segData';
% for home machine
wavfiledir='E:\Brandeis datasets\FMR1 Project Data\USV data\raw\Wav files';
outdir='E:\Brandeis datasets\FMR1 Project Data\USV data\segData';
wavfilenames=getAllFiles(wavfiledir);
fileinfo=dir(wavfiledir); fileinfo([fileinfo.isdir])=[];

if ~isfolder(outdir)
    mkdir(outdir);
end
nameadd='_callData';
bigclock=tic;
wb=waitbar(0,'Starting to build each session');
for i=1:length(wavfilenames)

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
    %{
    figure('Position',[2210 321 560 420]); subplot(2,1,1);
    imagesc(params.tvec(params.tvec>onoffsetm(1)-1 & params.tvec<onoffsetm(1,2)+1),...
        params.fvec, spect(:,params.tvec>onoffsetm(1)-1 & params.tvec<onoffsetm(1,2)+1));
    subplot(2,1,2);
    imagesc(params.tvec(params.tvec>onoffsetm(1)-1 & params.tvec<onoffsetm(1,2)+1),...
        params.fvec, thrshd(:,params.tvec>onoffsetm(1)-1 & params.tvec<onoffsetm(1,2)+1));
    yyaxis right;
    plot([onoffsetm(1) onoffsetm(1,2)],[.5 .5],'r','LineWidth',4)


    [sessdata.sessnotes]=input('Does this look okay (y/n)','s');
    close(gcf);
    %}
    segCalls=table(onoffset(:,1),onoffset(:,2),ones(size(onoffset,1),1),...
        'VariableNames',{'onsetTime','offsetTime','Accept'});
    
    callStats=getCallStats(segCalls,params,spect,thrshd);
    
    % reduce the data size of the spect
    spect=uint8(rescale(spect)*256);
    fprintf('%s contained %d calls \n', sessdata.name, height(segCalls))

    save(fullfile(outdir,[sessdata.name(1:end-4) nameadd]),...
        'sessdata','params','spect','blobs','segCalls','callStats');
    waitbar(i/length(wavfilenames),wb,sprintf('Its taken %d minutes, Likely %d more',...
        round(toc(bigclock)/60),round(toc(bigclock)/i*(length(wavfilenames)-i)/60)));
end
close(wb);
