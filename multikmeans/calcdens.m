%% Density data from signal and background
selected_signal = ones(longside,shortside);
to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');
[px, py] = gradient(image_MT_gray);
% Uniformity
Uniformity(k) = 100 * (1 - sqrt(var(to_analyse_all.PixelValues))/mean(to_analyse_all.PixelValues));
%Sparseness
Spars(k) = calcSparseness(to_analyse_all.PixelValues/mean(to_analyse_all.PixelValues),1);

% Kurtosis and skewness
skew(k) = skewness(to_analyse_all.PixelValues);

%Sdq - the root mean square gradient
Sdq(k) = sqrt(sum(px(:).*px(:)+(py(:).*py(:)))/length(to_analyse_all.PixelValues));
%Sdr - the developed interfacial area ratio

Sdr(k) = (sqrt(1+sum(px(:).*px(:)+(py(:).*py(:))))-1)/length(to_analyse_all.PixelValues);
% kmeans threshold

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
%bundling
mts_bundling(k) = (mean(image_MT_gray(signal_original~= 0),1)-...
    min(signal_original(signal_original>0))) / ...
    min(signal_original(signal_original>0));
%Spaces
im_bin_b = imcomplement(im_bin_c);
ccbg = bwconncomp(im_bin_b);
Space1 = regionprops(ccbg, 'Area');
if isempty(Space1)==0
    space(k) = mean(cat(1, Space1.Area));
else
    space(k) =0;
end

close all