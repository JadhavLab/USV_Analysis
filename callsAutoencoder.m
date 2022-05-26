

%%
% use create_tsne callback to generate this code

edit generate_VAE_encoder
edit generate_betaVAE_encoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% USE [encoderNet,decoderNet,data] = generate_VAE_encoder(imageStack)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now run the encoder decoder...
% Load the network model

% then yu can run tsne or umap on this and cluster
% first will need to pull a random set of usvs to generate our model.  This
% is throttled by the size of the gpu i guess...

% 2. build encoder and decoder using those images
%[encoderNet,decoderNet,data] = generate_VAE_encoder(allcalls);
% this needs work
%[encoderNet,decoderNet,data] = generate_betaVAE_encoder(allcalls);

%%
% make sure the vae works
nRecons=30;
images = dlarray(allcalls, 'SSCB');

  peekVAE(encoderNet,decoderNet,images,nRecons)
%

% 4. put the callduration back in
% 32 dims, 33rd is callduration, and allCallinfo is session and call number
% from runSess;
callData=[data zscore(allCalldurs(:,1))];

% 5. Show how these relate to age...
% callData=[32 hidden dims call length];

%
% two approaches, first get session averages and run that, second is run on
% all data across

figure;

for i=1:size(callData,2)
subplot(ceil(size(callData,2)),2,(i*2-1));
daymeans=accumarray(allCallinfo.sessionNumber,callData(:,i),[],@mean);
violindata={};
for k=1:max(runSess.age)
    violindata{k}=daymeans(runSess.age==k)';
end
agemeans=accumarray(runSess.age,daymeans,[],@mean);

agestd=accumarray(runSess.age,daymeans,[],@std);
okages=accumarray(allCallinfo.Age,1);

errorbar(find(okages>0),agemeans(okages>0),agestd(okages>0));
set(gca,'XLim',[3.5 20.5]);
subplot(ceil(size(callData,2)),2,(i*2));
violin(violindata(okages>0)); set(gca,'XTick',[1:10],'XTickLabel',2:2:20);
legend('off');
end

%%

% and now I pick some favorites and show what parameter they are spanning


% 1, 2, 3, 9, 12, 14, 18, 21, 23, 25, 26, 29, 31
% for each of these build say 5 images across the side from high to low,
% and then plot them on the side;
ndims=12;
freqdata=linspace(15,90,126);
patchcolors=lines(3);
for k=[1:12]
    figure;
    subplot(5,3,[1 2 4 5 7 8 10 11 13 14]);
    agemeans=accumarray(allCallinfo.Age,callData(:,k),[],@mean);
    daymeans=accumarray(allCallinfo.sessionNumber,callData(:,k),[],@mean);
    agestd=accumarray(runSess.age,daymeans,[],@std);
    okages=accumarray(allCallinfo.Age,1);
    highs=agemeans(okages>0)+agestd(okages>0);
    lows=flipud(agemeans(okages>0)-agestd(okages>0));
    %errorbar(find(okages>0),agemeans(okages>0),agestd(okages>0));
    patch([find(okages>0); flipud(find(okages>0))], [highs; lows],...
        patchcolors(1,:),'FaceAlpha',.5,'LineStyle','none');
    imagerands=repmat(randn(1,ndims)/2,5,1); imagerands(:,k)=[-2.5 -1 0 1 2.5];
    hold on; plot(find(okages>0),agemeans(okages>0));
    xlabel('age (PND)'); ylabel(sprintf('hyperparam %d score',k));
    for m=1:5
        randomNoise = dlarray(reshape(imagerands(m,:),[1,1,ndims]),'SSCB');
        generatedImage = predict(decoderNet, randomNoise);
        generatedImage = extractdata(generatedImage);
        subplot(5,3,m*3);
        imagesc(1:126,freqdata,SmoothMat2(zscore(real(generatedImage(2:end-1,2:end-1))),[10 10],3));
        set(gca,'YDir','normal'); xlabel('Time'); set(gca,'XTick',[],'CLim',[0 1]);
        ylabel('Freq (khz)');
    end
end
   
% the problem here is that the parameters are nonsensical.  I think I need
% to use either real calls as examples or i need to use beta vae
%% first look is pca

dayData=[];
for k=1:size(callData,2)
    dayData(:,k)=accumarray(allCallinfo.sessionNumber,callData(:,k),[],@mean);
end

figure;
% plot this out for m, f, fx, ctrl
groups=[USVlegend.age ((lower(USVlegend.Genotype)=='fx')*2 + double(lower(USVlegend.Genotype)=='wt'))];


allcolors=parula(length(unique(allCallinfo.Age)));
dubcolors=permute(repmat(allcolors,[1 1 3]),[2 3 1]);
[a,b,c]=pca(dayData);
% look at all calls regardless of session
gscatter(b(:,1),b(:,2),groups,dubcolors(:,:)','ox+',5,'off');
allkids=get(gca,'Children');
legend(flipud(allkids(end-2:end)),{'Het','wt','fx'});

figure;
[a,b,c]=pca(callData);
% look at all calls regardless of session
figure; gscatter(b(:,1),b(:,2),allCallinfo.Age);

figure; scatter3(b(:,1),b(:,2),b(:,3),2,allcolors(allCallinfo.Age,:),'filled');

% this suggests maybe we can capture the average session
sessPCs=[];
for i=1:3
    sessPCs(:,i)=accumarray(allCallinfo.sessionNumber,b(:,i),[],@mean);
end

figure; scatter3(sessPCs(:,1),sessPCs(:,2),sessPCs(:,3),5,allcolors(runSess.age,:),'filled');

%% run tsne on this

rng('default') % for fair comparison
Y = tsne(callData,'Perplexity',10,'Distance','cosine');
figure;
gscatter(Y(:,1),Y(:,2),allCallinfo.Age)
title('Cosine, color=age')



[a,b,c]=unique(allCallinfo.Genotype);
figure;
gscatter(Y(:,1),Y(:,2),c)
title('Cosine, color=Genotype')

%% or umap

if ~contains(path,'C:\Users\Jadhavlab\Documents\gitRepos\USV_Analysis\Umapfx\umap')
    addpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\USV_Analysis\Umapfx\umap'));
end
[embed2,umap2,clusterid] = run_umap(callData);
rmpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\USV_Analysis\Umapfx\umap'));



%% ica on these guys



addpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\USV_Analysis\PCA-ICA'));

Zfica = fastICA(callData',3);
b=Zfica';
% look at all calls regardless of session
figure; gscatter(b(:,1),b(:,2),allCallinfo.Age);

figure; scatter3(b(:,1),b(:,2),b(:,3),2,allcolors(allCallinfo.Age,:),'filled');

% this suggests maybe we can capture the average session
sessPCs=[];
for i=1:3
    sessPCs(:,i)=accumarray(allCallinfo.sessionNumber,b(:,i),[],@mean);
end

figure; scatter3(sessPCs(:,1),sessPCs(:,2),sessPCs(:,3),5,allcolors(runSess.age,:),'filled');


%% lets see if we can separate controls from fx

% first grab males:


for i=1:size(callData,2)
    subplot(size(callData,2),1,i);
    
    daymeans=accumarray(allCallinfo.sessionNumber,callData(:,i),[],@mean);
    cohortID=lower(runSess.Sex)=='f';
    for k=1:2
        
        agemeans=accumarray(runSess.age(cohortID==k-1),daymeans(cohortID==k-1),[],@mean);
        
        agestd=accumarray(runSess.age(cohortID==k-1),daymeans(cohortID==k-1),[],@std);
        okages=accumarray(runSess.age(cohortID==k-1),1);
        
        errorbar(find(okages>0),agemeans(okages>0),agestd(okages>0));
        hold on
        set(gca,'XLim',[3.5 20.5]);
    end
    legend('m','f');
end
% first grab males:
genos={'wt','het'};

for i=1:size(callData,2)
    subplot(size(callData,2),1,i);
    
    daymeans=accumarray(allCallinfo.sessionNumber,callData(:,i),[],@mean);
    okdays=lower(runSess.Sex)=='f';
    for k=1:2
        cohortID=(okdays & lower(runSess.Genotype)==genos{k});
        agemeans=accumarray(runSess.age(cohortID),daymeans(cohortID),[],@mean);
        
        agestd=accumarray(runSess.age(cohortID),daymeans(cohortID),[],@std);
        okages=accumarray(runSess.age(cohortID),1);
        
        errorbar(find(okages>0),agemeans(okages>0),agestd(okages>0));
        hold on
        set(gca,'XLim',[3.5 20.5]);
    end
end
