%% OVERLAY CALLS

% Uses ClusteringData cell from an Extracted Countours file
% Isolates all image maps for each call cluster
% Resizes and averages images to obtain an averaged call for each cluster
% Used to assess homeogeneity/validity of the clusters

resize = [300 600];
clusters = unique(T.Label);

figure(1);

for i=1:length(clusters)
    %Resize all the images that correspond to label i
    indexes = T.Label==i;
    rows_orig = find(indexes==1);
    rows = randsample(rows_orig,15); %Overlay a random subset of the images
    composite_im = uint8([zeros(300) zeros(300)]);
    for j=1:length(rows)
        subTCell = ClusteringData{rows(j),1};
        new_image = imresize(subTCell,resize);
        ClusteringData{rows(j),1}=new_image;
        composite_im = imfuse(composite_im,new_image);
    end
    subplot(4,6,i)
    image(composite_im);
    set(gca,'XTick',[], 'YTick', [])
    box off;
end