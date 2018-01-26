clear variables
clc
currdir = pwd;
addpath(pwd);
warning('off','stats:kmeans:FailedToConvergeRep');
cells = 10;
MTnumber = 25;
I = 25;
setup;

for Eccentricity = 0.7:0.1:0.9
    setcell;   
    for distribution = 30:10:40
        bundling = 0;
        for MTnumber = 25:25:150
            threshold = zeros(length(cells),1);
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
                
                clust2 = zeros((longside-2) * (shortside-2),14);
                parfor cl = 2:15
                    clust2(:,cl-1) = kmeans(Mat, cl, 'replicate',5);
                end
                eva = evalclusters(Mat,clust2,'DaviesBouldin');
                km = kmeans(Mat,eva.OptimalK,'replicate',5);
                km2 = reshape(km, longside-2,shortside-2);
                Image2_small = image_MT_gray(2:end-1,2:end-1);
                thr = zeros(eva.OptimalK,1);
                for clust = 1:eva.OptimalK
                    thr(clust) = mean(Image2_small(km2==clust));
                end
                [Num1, Idx1] = min(thr);
                threshold(k) = max(Image2_small(km2==Idx1));
                im_bin_c = imbinarize(image_MT_gray,threshold(k));
                
                calcdens;
                calcMTSD;
                
                
            end
            cd(currdir);
            %% Von Mises
            vonmises_fit_dist_sum;
            allcells;
        end
        MTnumber = 100;
        for bundling = 0:10:40    
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
                
                clust2 = zeros((longside-2) * (shortside-2),14);
                parfor cl = 2:15
                    clust2(:,cl-1) = kmeans(Mat, cl, 'replicate',5);
                end
                eva = evalclusters(Mat,clust2,'DaviesBouldin');
                km = kmeans(Mat,eva.OptimalK,'replicate',5);
                km2 = reshape(km, longside-2,shortside-2);
                Image2_small = image_MT_gray(2:end-1,2:end-1);
                thr = zeros(eva.OptimalK,1);
                for clust = 1:eva.OptimalK
                    thr(clust) = mean(Image2_small(km2==clust));
                end
                [Num1, Idx1] = min(thr);
                threshold(k) = max(Image2_small(km2==Idx1));
                im_bin_c = imbinarize(image_MT_gray,threshold(k));
                
                calcdens;
                calcMTSD;
                
                
            end
            cd(currdir);
            %% Von Mises
            vonmises_fit_dist_sum;
            allcells;
        end
    end
end
clear variables
close all
clc