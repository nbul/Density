%% Density data from signal and background
selected_signal = ones(longside,shortside);
to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');

% Skewness
skew(k) = skewness(to_analyse_all.PixelValues);

clust2 = zeros(longside * shortside,14);
ImageX = image_MT_gray(:);
parfor cl = 2:20
    clust2(:,cl-1) = kmeans(ImageX, cl, 'replicate',5);
end
eva = evalclusters(image_MT_gray(:),clust2,'DaviesBouldin');
km = kmeans(image_MT_gray(:),eva.OptimalK,'replicate',5);
Nculsters(k) = eva.OptimalK;
km2 = reshape(km, longside,shortside);
thr = zeros(eva.OptimalK,1);
for clust = 1:eva.OptimalK
    thr(clust) = mean(image_MT_gray(km==clust));
end
[Num1, Idx1] = min(thr);
threshold(k) = max(image_MT_gray(km2==Idx1));
im_bin_c = imbinarize(image_MT_gray,threshold(k));
signal_original = image_MT_gray .* im_bin_c;
% Signal Area
mts_area(k) = 100*sum(signal_original(:)>0)/longside/shortside;
close all