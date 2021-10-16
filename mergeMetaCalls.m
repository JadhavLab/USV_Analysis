%mergeMetaCalls

%% so first you need to load the metadata, and then you need to designate
% the call folder

load('C:\Users\Jadhavlab\Downloads\USVmetadata2021-10-15.mat');

callFolder='G:\USV data\Detections';
allCallFiles=cellfun(@(a) string(a), getAllFiles(callFolder,'.mat'));
%% attach the calls to each session
missingfiles=[];
for i=1:height(USVSession)
    fileName= USVSession.FileName(i);
    
    mymatch=contains(allCallFiles,append(fileName, " "));
    fprintf('File pattern is %s \n', fileName);
    fprintf('    Match File: %s \n',allCallFiles(mymatch));
   
    if sum(mymatch)~=1
        missingfiles=[missingfiles; fileName];
        fprintf( '   ***NO MATCH FILE FOR %s \n', fileName);
    else
        sessdat=load(allCallFiles(mymatch));
        USVSession.Calls{i}=sessdat.Calls;
        USVSession.audiodata{i}=sessdat.audiodata;
    end
    fprintf('\n');
end
        
    
    
    