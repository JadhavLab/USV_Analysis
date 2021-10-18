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
    
    mymatch=find(contains(allCallFiles,append(fileName, " ")));
    fprintf('File pattern is %s \n', fileName);
    fprintf('    Match File: %s \n',allCallFiles(mymatch));
   
    % try to fix
    if length(mymatch)~=1
        fprintf( '   ***NO MATCH FILE FOR %s \n', fileName);
        mymatch=input(sprintf('Match Filename in AllCallFiles with %s \n',...
            fileName));
    end
    % if fixed
    if length(mymatch)==1
        sessdat=load(allCallFiles(mymatch));
        USVSession.Calls{i}=sessdat.Calls;
        USVSession.audiodata{i}=sessdat.audiodata;
    else
        missingfiles=[missingfiles; fileName];
        
    end
    fprintf('\n');
end
        
%% attach stats to each session

statFolder='G:\USV data\CallStats2021-10-15';
allStatFiles=cellfun(@(a) string(a), getAllFiles(statFolder,'.xlsx'));
missingfiles=[];
for i=1:height(USVSession)
     fileName= USVSession.FileName(i);
    % add following space
     mymatch=find(contains(allStatFiles,append(fileName, " ")));
     fprintf('File pattern is %s \n', fileName);
     fprintf('    Match File: %s \n',allStatFiles(mymatch));
     % try to fix
    if length(mymatch)~=1
        fprintf( '   ***NO MATCH FILE FOR %s \n', fileName);
        mymatch=input(sprintf('Match Filename in AllCallFiles with %s \n',...
            fileName));
    end
    % if fixed
    if length(mymatch)==1
        callStats = importCallStats(allStatFiles(mymatch));
        USVSession.CallStats{i}=callStats;
    else
        missingfiles=[missingfiles; fileName];
    end
end