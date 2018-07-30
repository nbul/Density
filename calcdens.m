%% Density data from signal and background
selected_signal = ones(longside,shortside);
im_bin_c = imbinarize(image_MT_gray,mkthr(k));
im_bin_b = imcomplement(im_bin_c);

signal_original = image_MT_gray .* im_bin_c;
background_original = image_MT_gray .* (ones(longside,shortside) - im_bin_c);

to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');
to_analyse_o = regionprops(selected_signal, signal_original,'PixelValues');
to_analyse_c = regionprops(selected_signal, im_bin_c,'PixelValues');
to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');

sum_pixvalues_o = sum(to_analyse_o.PixelValues(:,1));
sum_pixvalues_back_o = sum(to_analyse_back_o.PixelValues(:,1));
num_pixvalues_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) ~= 0,1));
num_pixvalues_back_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
%Cell intensity
intensity(k) = mean(to_analyse_all.PixelValues);
%MTs intensity
intensity_mts(k) = mean(to_analyse_all.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
% Signal Area
mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);
% density
mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
    (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
%bundling
mts_bundling(k) = (mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))-...
    min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))) / ...
    min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
% Uniformity
Uniformity(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))./...
    (to_analyse_all.PixelValues + mean(to_analyse_all.PixelValues)))/length(to_analyse_all.PixelValues));
UNAAD(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))/...
    mean(to_analyse_all.PixelValues))/length(to_analyse_all.PixelValues));
%Sparseness
Spars(k) = calcSparseness(to_analyse_o.PixelValues/mean(to_analyse_o.PixelValues(to_analyse_o.PixelValues>0)),1);

% Kurtosis and skewness
if max(to_analyse_c.PixelValues)~= 0
    signal = to_analyse_all.PixelValues - min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
    kurt(k) = kurtosis(signal);
    skew(k) = skewness(signal);
else
    signal = 0;
    kurt(k) = 0;
    skew(k) = 0;
end

% Spaces
ccbg = bwconncomp(selected_signal.*im_bin_b);
Space1 = regionprops(ccbg, 'Area');
if isempty(Space1)==0
    space(k) = mean(cat(1, Space1.Area));
else
    space(k) =0;
end

Ent1(k) = entropy(double(mat2gray(image_MT_gray(selected_signal~= 0))));
Ent2(k) = entropy(double(mat2gray(image_MT_gray((selected_signal.*im_bin_c)~= 0))));

%Sdq - the root mean square gradient
[px, py] = gradient(image_MT_gray);
Sdq(k) = sqrt(sum(px(:).*px(:)+(py(:).*py(:)))/length(to_analyse_all.PixelValues));
%Sdr - the developed interfacial area ratio

Sdr(k) = (sqrt(1+sum(px(:).*px(:)+(py(:).*py(:))))-1)/length(to_analyse_all.PixelValues);

close all