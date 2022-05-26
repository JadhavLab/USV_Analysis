function peekVAE(encoderNet,decoderNet,XTest,nRecons)
%UNTITLED Summary of this function goes here



visualizeReconstruction(XTest,nRecons,encoderNet, decoderNet);

visualizeLatentSpace(XTest, encoderNet);

end




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

