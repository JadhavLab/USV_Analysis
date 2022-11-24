function [encoderNet,decoderNet,data] = generate_VAE_encoder(imageStack,encoderNet,decoderNet,verbose)
% Stolen somewhat from deepsqueak, but its also just using the matlab how
% to Variational AutoEncoder function

% data are
allcalls=imageStack;
numEpochs = 200;
miniBatchSize = 512;
nLatentDim=12;


% if we dont already have our net, build one from a random bit of these
% data

if ~exist('encoderNet','var') || ~exist('decoderNet','var')
    encoderNet=[]; decoderNet=[];
end

if ~exist('verbose','var'), verbose=0; end

images = dlarray(allcalls, 'SSCB');
    % Divide the images into training and validation
    [trainInd,valInd] = dividerand(size(images,4), .9, .1);
    XTrain  = images(:,:,:,trainInd); % 90% of data
    XTest   = images(:,:,:,valInd); % 10% of data

if isempty(encoderNet) ||isempty(decoderNet)
    % make sure calls are xpix by y pix by 1 by ncalls (third dim is color)
    [encoderNet, decoderNet] = VAE_model(nLatentDim);    
    % Train the network
    [encoderNet, decoderNet] = train_vae(encoderNet, decoderNet, XTrain, XTest,miniBatchSize,numEpochs);
end


if ~isempty(verbose) && verbose~=0
    % lets show the encoder and decoder reconstruct some images first
    visualizeReconstruction(XTest,10, encoderNet, decoderNet);
end


% now extract the low dimensional data
% [~, zMean] = sampling(encoderNet, images);
% zMean = stripdims(zMean)';
% zMean = gather(extractdata(zMean));
% data = double(zMean);

batchinds=[0 miniBatchSize:miniBatchSize:size(images,4)];
batches=[1+batchinds' [batchinds(2:end) size(images,4)]'];
data=[];
for bt=1:size(batches,1)
    [~,zMean]=sampling(encoderNet, images(:,:,:,batches(bt,1):batches(bt,2)));
    zMean = stripdims(zMean)';
    zMean = gather(extractdata(zMean));
    data = [data; double(zMean)];
end



end


function [encoderNet, decoderNet] = VAE_model(nLatentDim)

if ~exist('nLatentDim','var')
    latentDim = 32;
else
    latentDim=nLatentDim;
end
imageSize = [128, 128, 1];

encoderLG = layerGraph([
    imageInputLayer(imageSize,'Name','input_encoder','Normalization','none')
    
    convolution2dLayer(3, 8, 'Padding','same', 'Stride', 2, 'Name', 'conv1')
    batchNormalizationLayer('Name', 'bnorm1')
    reluLayer('Name','relu1')
    
    convolution2dLayer(3, 16, 'Padding','same', 'Stride', 2, 'Name', 'conv2')
    batchNormalizationLayer('Name', 'bnorm2')
    reluLayer('Name','relu2')
    
    convolution2dLayer(3, 32, 'Padding','same', 'Stride', 2, 'Name', 'conv3')
    batchNormalizationLayer('Name', 'bnorm3')
    reluLayer('Name','relu3')
    
    convolution2dLayer(3, 64, 'Padding','same', 'Stride', 2, 'Name', 'conv4')
    batchNormalizationLayer('Name', 'bnorm4')
    reluLayer('Name','relu4')
    
    fullyConnectedLayer(1024, 'Name', 'fc_1')
    reluLayer('Name','relu5')

    fullyConnectedLayer(2 * latentDim, 'Name', 'fc_encoder')
    ]);

decoderLG = layerGraph([
    imageInputLayer([1 1 latentDim],'Name','i','Normalization','none')
    
    transposedConv2dLayer(16, 32, 'Cropping', 0, 'Stride', 1, 'Name', 'transpose1')
    batchNormalizationLayer('Name', 'bnorm1')
    reluLayer('Name','relu1')
    transposedConv2dLayer(3, 32, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose2')
    batchNormalizationLayer('Name', 'bnorm2')
    reluLayer('Name','relu2')
    transposedConv2dLayer(3, 24, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose3')
    batchNormalizationLayer('Name', 'bnorm3')
    reluLayer('Name','relu3')
    transposedConv2dLayer(3, 16, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose4')
    batchNormalizationLayer('Name', 'bnorm4')
    reluLayer('Name','relu4')
    transposedConv2dLayer(3, 8, 'Cropping', 'same', 'Stride', 1, 'Name', 'transpose5')
    batchNormalizationLayer('Name', 'bnorm5')
    reluLayer('Name','relu5')
    transposedConv2dLayer(3, 1, 'Cropping', 'same', 'Name', 'transpose6')
    ]);



% analyzeNetwork(encoderLG)
% analyzeNetwork(decoderLG)

encoderNet = dlnetwork(encoderLG);
decoderNet = dlnetwork(decoderLG);

end


function [encoderNet, decoderNet] = train_vae(encoderNet, decoderNet, XTrain, XTest,miniBatchSize,numEpochs)

numTrainImages = size(XTrain, 4);

executionEnvironment = "auto";

lr = 1e-3;
numIterations = floor(numTrainImages/miniBatchSize);
iteration = 0;

avgGradientsEncoder = [];
avgGradientsSquaredEncoder = [];
avgGradientsDecoder = [];
avgGradientsSquaredDecoder = [];


figure1 = figure('Color',[1 1 1],'Position',[200 200 600 500]);
axes1 = axes('Parent',figure1,'LineWidth',1,'TickDir','out',...
    'FontSmoothing','on',...
    'FontSize',12);
ylabel(axes1,'ELBO loss');
xlabel(axes1,'Epoch');
plotTitle = title(axes1, 'Training progress', 'Close this window to end training');
h = animatedline(axes1, 'Color', [.1, .9, .7], 'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 20);
% xlim(axes1, [0, numEpochs])

set(axes1, 'yscale', 'log')
for epoch = 1:numEpochs
    tic;
    % this runs 200 (numIterations) iterations of network updating to
    % generate a new set of y (x=images, y=sparse vector)
    for i = 1:numIterations
        iteration = iteration + 1;
        idx = (i-1)*miniBatchSize+1:i*miniBatchSize;
        XBatch = XTrain(:,:,:,idx);
        XBatch = dlarray(single(XBatch), 'SSCB');
        
        if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
            XBatch = gpuArray(XBatch);           
        end 
            
        [infGrad, genGrad] = dlfeval(...
            @modelGradients, encoderNet, decoderNet, XBatch);
        
        [decoderNet.Learnables, avgGradientsDecoder, avgGradientsSquaredDecoder] = ...
            adamupdate(decoderNet.Learnables, ...
                genGrad, avgGradientsDecoder, avgGradientsSquaredDecoder, iteration, lr);
        [encoderNet.Learnables, avgGradientsEncoder, avgGradientsSquaredEncoder] = ...
            adamupdate(encoderNet.Learnables, ...
                infGrad, avgGradientsEncoder, avgGradientsSquaredEncoder, iteration, lr);
    end
    elapsedTime = toc;
    
    % if there is no figure, return
    if ~isvalid(h)
        return
    end

    
    
    % forward the whole set
    [z, zMean, zLogvar] = sampling(encoderNet, XTest);
    forward(encoderNet, XTest);
    % backward, generate the predictions
    % and this is where we run into problems
    % lets grab this in sections
    batchinds=[0 miniBatchSize:miniBatchSize:size(z,4)];
    batches=[1+batchinds' [batchinds(2:end) size(z,4)]'];
    xPred=[];
    for bt=1:size(batches,1)
        xPredRaw=forward(decoderNet, z(:,:,:,batches(bt,1):batches(bt,2)));
        xPred = cat(4,xPred,sigmoid(xPredRaw));
    end
        
   
    elbo = ELBOloss(XTest, xPred, zMean, zLogvar);
    
    % Update figure and print results
    fprintf('Epoch : %-3g Test ELBO loss = %#.5g. Time taken for epoch = %#.3gs\n', epoch, gather(extractdata(elbo))/2, elapsedTime)
    addpoints(h,epoch,double(gather(extractdata(elbo))));
    plotTitle.String = sprintf('Training progress - epoch %u/%u', epoch, numEpochs);
    drawnow 
end


end

function [infGrad, genGrad] = modelGradients(encoderNet, decoderNet, x)

% generate your z distribution to get your KL values
[z, zMean, zLogvar] = sampling(encoderNet, x);
% generate your Xpred for reconstruction loss
xPred = sigmoid(forward(decoderNet, z));
% merge the two
loss = ELBOloss(x, xPred, zMean, zLogvar);
% and run gradient
[genGrad, infGrad] = dlgradient(loss, decoderNet.Learnables, ...
    encoderNet.Learnables);

end


% this generates your z distributions
function [zSampled, zMean, zLogvar] = sampling(encoderNet, x)
compressed = forward(encoderNet, x);
d = size(compressed,1)/2;
zMean = compressed(1:d,:);
zLogvar = compressed(1+d:end,:);

sz = size(zMean);
epsilon = randn(sz); % get rand normally distributed (will be standard z)
sigma = exp(.5 * zLogvar); % get your variance
z = epsilon .* sigma + zMean; % basically reshape your z distrib
z = reshape(z, [1,1,sz]); % reshape these variables by dimension
% this basically will allow you to generate a sample of data points, so
% that you can get a conditional distribution of real values, given a
% distribution of answers, this is the bayes part

zSampled = dlarray(z, 'SSCB'); % send into your dlarray
end



%%%%% LOSS FUNCTION STUFF %%%%%

% elbo loss: −LVAE=logpθ(x)−DKL(qϕ(z|x)∥pθ(z|x))

% when b=1 its the normal ELBO function. when b>1, it splits the bariables.
%
% beta vae loss: LBETA(ϕ,β)=−Ez∼qϕ(z|x)logpθ(x|z)+βDKL(qϕ(z|x)∥pθ(z))
% where the addition is the Ez~qw(z|x)


% this is the loss function, 
% which is where you would compute beta

% this estimates the kullback-leibler divergence
%  math: Dkl= real distrib (imges given vector) ||(similarity) estimated
%  distrib(images given vector)
% the key here is that you can calculate the real



function elbo = ELBOloss(x, xPred, zMean, zLogvar)

% xpred=images, zMean=latent variables, x is real images and zlogvar is how
% you generate your posterior distributions (gaussian with var Z)
squares = 0.5*(xPred-x).^2;
reconstructionLoss  = sum(squares, [1,2,3]); % e.g. likelihood of generating the data you did


KL = -.5 * sum(1 + zLogvar - zMean.^2 - exp(zLogvar), 1); % kl which is 

% elbo= evidence lower bound (average across your test image batch)
elbo = mean(reconstructionLoss + KL);
end






function visualizeReconstruction(XTest,nRecons, encoderNet, decoderNet)

for c=1:nRecons
    idx = randi(size(XTest,4),1); % pull random
    X = XTest(:,:,:,idx);
    
    [z, ~, ~] = sampling(encoderNet, X);
    XPred = sigmoid(forward(decoderNet, z));
    
    X = gather(extractdata(X));
    XPred = gather(extractdata(XPred));
    
    comparison = [X, ones(size(X,1),1), XPred];
    figure; imshow(comparison,[]), title("Example ground truth image vs. reconstructed image")

end

end



function visualizeLatentSpace(XTest, encoderNet)
[~, zMean, zLogvar] = sampling(encoderNet, XTest);

zMean = stripdims(zMean)';
zMean = gather(extractdata(zMean));

zLogvar = stripdims(zLogvar)';
zLogvar = gather(extractdata(zLogvar));

[~,scoreMean] = pca(zMean);
[~,scoreLogvar] = pca(zLogvar);

c = parula(10);
f1 = figure;
figure(f1)
title("Latent space")

ah = subplot(1,2,1);
scatter(scoreMean(:,2),scoreMean(:,1),[]);
ah.YDir = 'reverse';
axis equal
xlabel("Z_m_u(2)")
ylabel("Z_m_u(1)")

ah = subplot(1,2,2);
scatter(scoreLogvar(:,2),scoreLogvar(:,1),[]);
ah.YDir = 'reverse';
xlabel("Z_v_a_r(2)")
ylabel("Z_v_a_r(1)")
axis equal
end

% inputs- decoderNet, and latentDim is number of latent dimensions, and
% this produces 25 images
function generate(decoderNet, latentDim)
randomNoise = dlarray(randn(1,1,latentDim,1),'SSCB');
generatedImage = sigmoid(predict(decoderNet, randomNoise));
generatedImage = extractdata(generatedImage);

f3 = figure;
figure(f3)
imshow(imtile(generatedImage, "ThumbnailSize", [100,100]))
title("Generated random samples")
drawnow
end

%{
notes
this is the scripting that was used to generate this function.


[encoderNet, decoderNet] = VAE_model();

% data are

% build encoder from a set stack size... not sure how much it can handle...
randpull=randperm(size(allCallinfo,1));
pullsize=20000; % start with twenty thousand
imageInfo=allCallinfo(randpull(1:pullsize),:);


images = dlarray(allcalls(:,:,:,randpull(1:pullsize)), 'SSCB');



% Divide the images into training and validation
[trainInd,valInd] = dividerand(size(images,4), .9, .1);
XTrain  = images(:,:,:,trainInd);
XTest   = images(:,:,:,valInd);

% Train the network
[encoderNet, decoderNet] = train_vae(encoderNet, decoderNet, XTrain, XTest);

% now pull all the image data
myblocks=[[1 pullsize+1:pullsize:size(allCallinfo,1)]' [pullsize:pullsize:size(allCallinfo,1) size(allCallinfo,1)]' ];
alldata=[];
for i=1:size(myblocks,1)
% now extract the low dimensional data
[~, zMean] = sampling(encoderNet, images(:,:,:,myblocks(i,1):myblocks(i,2)));
zMean = stripdims(zMean)';
zMean = gather(extractdata(zMean));
data = double(zMean);
alldata=[alldata; data];
end
% and we already have the original imageinfo
allimageShort=[alldata imageInfo(:,1)];

% now add classes (age, or genotype) and image it




% now reform this and get some preliminary results?
data=[data imageInfo(:,1)]; % tack on call duration
%}