function [encoderNet,decoderNet,data] = generate_VAE_encoder(imageStack)
% Stolen somewhat from deepsqueak, but its also just using the matlab how
% to

% make sure calls are xpix by y pix by 1 by ncalls (third dim is color)
[encoderNet, decoderNet] = VAE_model();

% data are
allcalls=imageStack;

% take random values
randpull=randperm(size(allcalls,4));
pullsize=10000; % start with ten thousand
%imageInfo=allCallinfo(randpull(1:pullsize),:);


images = dlarray(allcalls(:,:,:,randpull(1:pullsize)), 'SSCB');

% Divide the images into training and validation
[trainInd,valInd] = dividerand(size(images,4), .9, .1);
XTrain  = images(:,:,:,trainInd);
XTest   = images(:,:,:,valInd);

% Train the network
[encoderNet, decoderNet] = train_vae(encoderNet, decoderNet, XTrain, XTest);

% now pull 

% now extract the low dimensional data
[~, zMean] = sampling(encoderNet, images);
zMean = stripdims(zMean)';
zMean = gather(extractdata(zMean));
data = double(zMean);



end


function [encoderNet, decoderNet] = VAE_model()


latentDim = 32;
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


function [encoderNet, decoderNet] = train_vae(encoderNet, decoderNet, XTrain, XTest)

numTrainImages = size(XTrain, 4);

executionEnvironment = "auto";


numEpochs = 200;
miniBatchSize = 128;
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
    
    if ~isvalid(h)
        return
    end
    
    [z, zMean, zLogvar] = sampling(encoderNet, XTest);
    forward(encoderNet, XTest);
    xPred = sigmoid(forward(decoderNet, z));
    elbo = ELBOloss(XTest, xPred, zMean, zLogvar);
    
    % Update figure and print results
    fprintf('Epoch : %-3g Test ELBO loss = %#.5g. Time taken for epoch = %#.3gs\n', epoch, gather(extractdata(elbo))/2, elapsedTime)
    addpoints(h,epoch,double(gather(extractdata(elbo))));
    plotTitle.String = sprintf('Training progress - epoch %u/%u', epoch, numEpochs);
    drawnow 
end
end

function [zSampled, zMean, zLogvar] = sampling(encoderNet, x)
compressed = forward(encoderNet, x);
d = size(compressed,1)/2;
zMean = compressed(1:d,:);
zLogvar = compressed(1+d:end,:);

sz = size(zMean);
epsilon = randn(sz);
sigma = exp(.5 * zLogvar);
z = epsilon .* sigma + zMean;
z = reshape(z, [1,1,sz]);
zSampled = dlarray(z, 'SSCB');
end

function [infGrad, genGrad] = modelGradients(encoderNet, decoderNet, x)
[z, zMean, zLogvar] = sampling(encoderNet, x);
xPred = sigmoid(forward(decoderNet, z));
loss = ELBOloss(x, xPred, zMean, zLogvar);
[genGrad, infGrad] = dlgradient(loss, decoderNet.Learnables, ...
    encoderNet.Learnables);
end

function elbo = ELBOloss(x, xPred, zMean, zLogvar)
squares = 0.5*(xPred-x).^2;
reconstructionLoss  = sum(squares, [1,2,3]);

KL = -.5 * sum(1 + zLogvar - zMean.^2 - exp(zLogvar), 1);

elbo = mean(reconstructionLoss + KL);
end

function visualizeReconstruction(XTest,nRecons, encoderNet, decoderNet)

title("Example ground truth image vs. reconstructed image")
for c=1:nRecons
    idx = randi(size(XTest,4),1); % pull random
    X = XTest(:,:,:,idx);
    
    [z, ~, ~] = sampling(encoderNet, X);
    XPred = sigmoid(forward(decoderNet, z));
    
    X = gather(extractdata(X));
    XPred = gather(extractdata(XPred));
    
    comparison = [X, ones(size(X,1),1), XPred];
    figure; imshow(comparison,[]),
end

end



function visualizeLatentSpace(XTest, YTest, encoderNet)
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
cb = colorbar; cb.Ticks = 0:(1/9):1; cb.TickLabels = string(0:9);

ah = subplot(1,2,2);
scatter(scoreLogvar(:,2),scoreLogvar(:,1),[]);
ah.YDir = 'reverse';
xlabel("Z_v_a_r(2)")
ylabel("Z_v_a_r(1)")
cb = colorbar;  cb.Ticks = 0:(1/9):1; cb.TickLabels = string(0:9);
axis equal
end

function generate(decoderNet, latentDim)
randomNoise = dlarray(randn(1,1,latentDim,25),'SSCB');
generatedImage = sigmoid(predict(decoderNet, randomNoise));
generatedImage = extractdata(generatedImage);

f3 = figure;
figure(f3)
imshow(imtile(generatedImage, "ThumbnailSize", [100,100]))
title("Generated samples of digits")
drawnow
end