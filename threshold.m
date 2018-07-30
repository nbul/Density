%% kmeans threshold
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
mkthr(k) = max(image_MT_gray(km2==Idx1));
levelsmk(k) = eva.OptimalK;
%% Otsu threshold
otsuthr(k) = graythresh(image_MT_gray/255)*255;

%% Edges threshold
image_edges = edge(image_MT_gray, 'Canny');
signal_edges = image_MT_gray .* image_edges;
[A1, A2] = histcounts(signal_edges(signal_edges>0));
bins = (min(A2)+((A2(2)-A2(1))/2)):((A2(2)-A2(1))):(max(A2)-((A2(2)-A2(1))/2));
[XOut,YOut] = prepareCurveData(bins,A1);
fo = fitoptions('gauss1', 'Lower', [0 min(A2) 0], 'Upper', [Inf max(A2) Inf]);
[thr, gof_edges] = fit(XOut, YOut, 'gauss1', fo);
thredges(k) = thr.b1;

%% threshold 2kmeans
image2 = image_MT_gray(:);
km = kmeans(image2,2,'replicate',5);
thr = zeros(2,1);
for clust = 1:2
    thr(clust) = mean(image2(km==clust));
end
[Num2, Idx2] = min(thr);
thr2kmeans(k) = max(image2(km==Idx2));

%% Multiotsu threshold
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
