%% OVERLAY CALLS

% Uses ClusteringData cell from an Extracted Countours file
% Isolates all image maps for each call cluster
% Resizes and averages images to obtain an averaged call for each cluster
% Used to assess homeogeneity/validity of the clusters

resize = [300 600];


indexes = ClusteringData{:,1}==;
subT = T(indexes,:);


new_images = imresize(ClusteringData{:,1},resize);