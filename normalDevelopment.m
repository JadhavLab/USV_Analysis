% USVlegend is saved in USVmetadataC1-9.at
% lets take just the males, and just between P4 to P10
% this is hardcoded for now, but suffice to say its going to be 600 msec
% images from 15 khz to 90 khz


runSess=USVlegend(lower(USVlegend.Sex)=='m' & USVlegend.age>13,:);
%runSess=USVlegend(USVlegend.age>10,:);
allcalls=[]; allCalldurs=[];
wb=waitbar(0,'concatenating animal images');
% what percentage of images will we need to classify these guys????
% probably like half
bigclock=tic;
for i=1:height(runSess)

    load(fullfile(runSess.folder(i,:),runSess.name(i)),'blobs','segCalls','params','spect');
    okfreqs=params.fvec>15000 & params.fvec<90000;
    badcalls=segCalls.onsetTime<5 | segCalls.offsetTime>min([180 size(blobs,2)*params.timestep])-5;
    segCalls=segCalls(~badcalls,:);
    mycalls=false(128,128,1,height(segCalls));
    callcenters=mean(table2array(segCalls(:,1:2)),2);
    centerinds=interp1(params.tvec,1:length(params.tvec),callcenters,'nearest');
    calldurs=(segCalls.offsetTime-segCalls.onsetTime)/params.timestep;
    

    for j=1:length(callcenters)
        % go forward 200 msec, and forward 200 msec, then resize to 256
        % were going from many to 128x128 x 1 images- this is the bw
        %mycalls(:,:,:,j)=imresize(blobs(okfreqs,(centerinds(j)-round(calldurs(j)/1.9)):(centerinds(j)+round(calldurs(j)/1.9))),[128,128]);
        
        % or use the actual image, which may be better...
        mycalls(:,:,:,j)=imresize(spect(okfreqs,(centerinds(j)-round(calldurs(j)/1.9)):(centerinds(j)+round(calldurs(j)/1.9))),[128,128]);

        % or use both! but you'll have to edit here!
        %mycalls(:,:,1,j)=imresize(blobs(okfreqs,(centerinds(j)-round(calldurs(j)/1.9)):(centerinds(j)+round(calldurs(j)/1.9))),[128,128]);
        %mycalls(:,:,2,j)=imresize(spect(okfreqs,(centerinds(j)-round(calldurs(j)/1.9)):(centerinds(j)+round(calldurs(j)/1.9))),[128,128]);

    end
    allCalldurs=[allCalldurs; [calldurs ones(size(calldurs,1),1)*i (1:length(calldurs))']]; % add call duration, then the sessnum and callnum
    allcalls=cat(4,allcalls,mycalls);
    waitbar(i/height(runSess),wb,sprintf('Has taken %d mins, will prolly go %d more mins',...
    round(toc(bigclock)/60),round(toc(bigclock)/i*(height(runSess)-i)/60)));
end
close(wb);


allCallinfo=table(allCalldurs(:,2),allCalldurs(:,3),'VariableNames',...
    {'sessionNumber','call number'});

for i=1:height(allCallinfo)
    allCallinfo.Age(i)=runSess.age(allCallinfo.sessionNumber(i));
    allCallinfo.Sex(i)=runSess.Sex(allCallinfo.sessionNumber(i));
    allCallinfo.Cohort(i)=runSess.Cohort(allCallinfo.sessionNumber(i));
    allCallinfo.Rat(i)=runSess.ratname(allCallinfo.sessionNumber(i));
    allCallinfo.Genotype(i)=runSess.Genotype(allCallinfo.sessionNumber(i));
end

%save('AllCallImages','runSess','allCallinfo','allcalls','-v7.3')

%%

% this is all images for all sessions:
load('G:\USV data\AllCallImages');


[data] = run_VAE_encoder(encoderNet,allcalls);
%%
% use create_tsne callback to generate this code

%edit generate_VAE_encoder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% USE [encoderNet,decoderNet,data] = generate_VAE_encoder(imageStack)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now run the encoder decoder...
% Load the network model

% then yu can run tsne or umap on this and cluster
% first will need to pull a random set of usvs to generate our model.  This
% is throttled by the size of the gpu i guess...


edit callsAutoencoder.m
