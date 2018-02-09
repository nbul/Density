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

%threshold otsu
throtsu(k) = graythresh(image_MT_gray/255)*255;
%threshold edges
image_edges = edge(image_MT_gray, 'Canny');
signal_edges = image_MT_gray .* image_edges;
[A1, A2] = histcounts(signal_edges(signal_edges>0));
bins = (min(A2)+((A2(2)-A2(1))/2)):((A2(2)-A2(1))):(max(A2)-((A2(2)-A2(1))/2));
[XOut,YOut] = prepareCurveData(bins,A1);
fo = fitoptions('gauss1', 'Lower', [0 min(A2) 0], 'Upper', [Inf max(A2) Inf]);
[thr, gof_edges] = fit(XOut, YOut, 'gauss1', fo);
thredges(k) = thr.b1;
%threshold 2kmeans
image2 = image_MT_gray(:);
km = kmeans(image2,2,'replicate',5);
thr = zeros(2,1);
for clust = 1:2
    thr(clust) = mean(image2(km==clust));
end
[Num1, Idx1] = min(thr);
thr2kmeans(k) = max(image2(km==Idx1));
%multithresh
thrtemp = zeros(19,1);
metric = zeros(19,1);
for T = 2:20
    [TT, TM] = multithresh(image_MT_gray, T );
    thrtemp(T-1) = min(TT);
    metric(T-1) = TM;
end
[NumT, IdxT] = max(metric);
thrmulti(k) = thrtemp(IdxT);
levelsmulti(k) = IdxT + 1;
close all