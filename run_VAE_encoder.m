function [data] = run_VAE_encoder(encoderNet,imageStack)




numEpochs = 200;
miniBatchSize = 512;
nLatentDim=6;


images = dlarray(imageStack, 'SSCB');

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

