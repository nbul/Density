clear variables
clc
currdir = pwd;
addpath(pwd);
warning('off','stats:kmeans:FailedToConvergeRep');
setup;
setcell;
for k=1:cells
    MT;
    % Threshold
    Mat = zeros((shortside-2)*(longside-2),3);
    im_adjusted = imadjust(image_MT_gray/255);
    counter4=0;
    for xc=2:(shortside-1)
        for yc=2:(longside-1)
            counter4=counter4+1;
            Mat(counter4,1) = image_MT_gray(yc,xc);
            Mat(counter4,2) = (image_MT_gray(yc-1,xc-1) + image_MT_gray(yc-1,xc) + image_MT_gray(yc-1,xc+1) +...
                image_MT_gray(yc+1,xc-1) + image_MT_gray(yc+1,xc) + image_MT_gray(yc+1,xc+1) +...
                image_MT_gray(yc,xc-1) + image_MT_gray(yc,xc+1))/8;
            median_mat = image_MT_gray(yc-1:yc+1,xc-1:xc+1);
            Mat(counter4,3) = median(median_mat(:));
        end
    end
    
    myfunc = @(X,K)(kmeans(X, K, 'replicate',5));
    eva = evalclusters(Mat,myfunc,'DaviesBouldin',...
        'klist',2:8);
    
    km = kmeans(Mat,eva.OptimalK,'replicate',5);
    km2 = reshape(km, longside-2,shortside-2);
    thr = zeros(eva.OptimalK,1);
    for clust = 1:eva.OptimalK
        thr(clust) = mean(image_MT(km2==clust));
    end
    [Num1, Idx1] = min(thr);
    threshold(k) = max(image_MT(km2==Idx1));
    im_bin_c = imbinarize(image_MT_gray,threshold(k));
    
    calcdens;
    calcMTSD;
    
 
end
cd(currdir);
%% Von Mises
vonmises_fit_dist_sum;
allcells;
clear variables
close all
clc