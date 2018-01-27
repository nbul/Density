bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
m_added_norm(:,1) = bincenter;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';
side = 50;
for angle = 1:89
    X = side/2;
    Y = side/2;
    a = 1/tand(angle);
    b = Y - a*X;
    intersect=zeros(1,4);
    l=0;
    X_temp = 1;
    Y_temp = ceil(b+a);
    if Y_temp>0 && Y_temp<=side
        l=1;
        intersect(1,1) = X_temp;
        intersect(1,2) = Y_temp;
    end
    
    X_temp = ceil(1 - b/a);
    Y_temp = 1;
    if X_temp>0 && X_temp<=side
        if l == 0
            l = 1;
            intersect(1,1) = X_temp;
            intersect(1,2) = Y_temp;
        else
            l = l+1;
            intersect(1,3) = X_temp;
            intersect(1,4) = Y_temp;
        end
    end
    
    X_temp = ceil((side - b)/a);
    Y_temp = side;
    if X_temp>0 && X_temp<=side
        if l == 0
            l = 1;
            intersect(1,1) = X_temp;
            intersect(1,2) = Y_temp;
        else
            l = l+1;
            intersect(1,3) = X_temp;
            intersect(1,4) = Y_temp;
        end
    end
    
    X_temp = side;
    Y_temp = a*side +b;
    if Y_temp>0 && Y_temp<=side
        if l == 0
            l = 1;
            intersect(1,1) = X_temp;
            intersect(1,2) = Y_temp;
        else
            l = l+1;
            intersect(1,3) = X_temp;
            intersect(1,4) = Y_temp;
        end
    end
    
    image = zeros(side,side);
    image_MT_gray = image;
    
    image_MT = insertShape(image,'line',intersect(1,:), 'linewidth', 3, 'Color', [50 50 50]);
    image_MT_gray = image_MT_gray + image_MT(:,:,1);
    
    for gaus = 1:4
        image_MT_gray = image_MT_gray + 10;
        image_MT_gray(image_MT_gray>255) = 255;
        if gaus>1
            image_MT_gray = imgaussfilt(image_MT_gray,gaus-1);
        end
        clear H_full V_full H V M D x y mxd_thr mxd_corrected mxd_indexed
        V = fig_to_LCN(image_MT_gray,3);
        
        figure(2)
            subplot(s1,s2,1); imshow(image_MT_gray, []); 
            subplot(s1,s2,2); imshow(V, []); 
        pause 

        object_double = V;
        H_full = conv2(object_double,Gx);
        V_full = conv2(object_double,Gy);
        H = H_full(5:side,5:side);
        V = V_full(5:side,5:side);
        M = sqrt(H.^2 + V.^2);
        D = -(180/pi) * atan2(V, H);
        
        [x, y] = size(M);
        p = 1;
        mxd = zeros(1,2);
        
        for j = 2:(y-1)
            for i=2:(x-1)
                %Only directions different to zero are added to table
                if ((M(i,j)) & (M(i+1,j)) & (M(i-1,j)) & (M(i,j+1))  & (M(i,j-1) ~=0) &...
                        (M(i-1,j-1)) & (M(i-1 , j+1)) & (M(i+1 , j-1)) & (M(i+1, j+1))) ~= 0
                    mxd(p,2) = M(i,j); %Second column with magnitudes
                    mxd(p,1) = D(i,j); %First column with angles
                    p = p + 1;
                end
            end
        end
        max_mxd  = max(mxd(:,2)); %maximum magnitude
        mxd_thr = mxd./repmat([1,max_mxd], length(mxd), 1); %normalised to max magnitude
        
        % Remove all pixels with magnitude less than 22% of maximum
        
        mxd_corrected = mxd_thr(mxd_thr(:,2) >= 0.22,:);
        mxd_corrected(mxd_corrected(:,1) < 0,1) = mxd_corrected(mxd_corrected(:,1) < 0,1) + 180;
        mxd_corrected(:,1) = mxd_corrected(:,1) - 90;
        mxd_corrected(mxd_corrected(:,1) >= 90,1) = 89.9;
        mxd_corrected = sortrows(mxd_corrected,1);
        
        % Make histogram
        [N, bins] = histc(mxd_corrected(:,1),binrange);
        
        mxd_indexed(:,2) = mxd_corrected(:,2);
        mxd_indexed(:,1) = bins(:,1);
        % Make distribution
        m_added = zeros(45,1);
        for i=1:(length(binrange)-1)
            m_added(i,1) = sum(mxd_indexed(mxd_indexed(:,1) == i,2));
        end
        
        m_added_norm(:,2) = m_added/sum(m_added);
        vonmises_fit_dist_sum;
        SD2(angle,gaus) = SD;
        mu2(angle,gaus) = mu+90;
    end
end
s1 = 1;
s2 = 2;
image1 = figure;
subplot(s1,s2,1);plot(2:89,mu2(2:89,1), 2:89, mu2(2:89,2), 2:89, mu2(2:89,3), 2:89, mu2(2:89,4))
subplot(s1,s2,2);plot(2:89,SD2(2:89,1), 2:89, SD2(2:89,2), 2:89, SD2(2:89,3), 2:89, SD2(2:89,4))

figure(2)
subplot(s1,s2,1); imshow(image_MT_gray, []); 
subplot(s1,s2,2); imshow(V, []);