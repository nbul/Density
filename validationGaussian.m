clear variables
clc

MTnumber = 100;
Eccentricity = 0.9;
shortside = 200;


distribution = 30;
I = 25;
bundling = 1;
method = 1;
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';

currdir = pwd;
addpath(pwd);
m_added_norm = zeros(45,11);
m_added_norm(:,1) = bincenter;
Value = zeros(45,11);
Value(:,1) = bincenter;
mts_density = zeros(1, 10);
mts_area = zeros(1, 10);
Uniformity = zeros(1, 10);
Spars = zeros(1, 10);
mts_bundling = zeros(1, 10);
kurt = zeros(1, 10);
skew = zeros(1, 10);

long_side = zeros(1, (4*shortside));
for i=(2*shortside):(10*shortside)
    selected_signal = ones(i,shortside);
    to_analyse_o = regionprops(selected_signal, selected_signal,'Eccentricity');
    long_side(i-2*shortside+1) = abs(to_analyse_o.Eccentricity-Eccentricity);
end
[Ecc, Ind] = min(long_side);
longside = Ind + 2*shortside-1;

Length = zeros(1,MTnumber);
result_dir = '/Users/nataliabulgakova/MT-project/Robustness/Densityvalidation';
cd('/Users/nataliabulgakova/MT-project/Robustness/Densityvalidation');
% image_dir_name = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
%     num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_', num2str(method)];
% if exist([result_dir,'/', image_dir_name],'dir') == 0
%     mkdir(result_dir,image_dir_name);
% end
% image_dir = [result_dir,'/', image_dir_name];
% cd(image_dir);
for G=0:4
    for k=1:10
        %% Generate random dots within the cell
        rng('shuffle');
        X = randi(shortside, MTnumber,1);
        Y = randi(longside, MTnumber,1);
        %% Generate angles with a given distribution
        angles = normrnd(0,distribution, [MTnumber 1]);
        for i=1:MTnumber
            if angles(i, 1) <= -90
                angles(i, 1) = angles(i, 1) + 180;
            elseif angles(i, 1) > 90
                angles(i, 1) = angles(i, 1) - 180;
            end
        end
        
        bundled = randi(bundling, MTnumber,1);
        %% Line parameters and start/end points
        a = 1./tand(angles);
        b = Y - a.*X;
        intersect=zeros(MTnumber,4);
        l=zeros(MTnumber,1);
        for i=1:MTnumber
            l(i)=0;
            X_temp = 1;
            Y_temp = ceil(b(i)+a(i));
            if Y_temp>0 && Y_temp<=longside
                l(i)=1;
                intersect(i,1) = X_temp;
                intersect(i,2) = Y_temp;
            end
            
            X_temp = ceil(1 - b(i)/a(i));
            Y_temp = 1;
            if X_temp>0 && X_temp<=shortside
                if l(i) == 0
                    l(i) = 1;
                    intersect(i,1) = X_temp;
                    intersect(i,2) = Y_temp;
                else
                    l(i) = l(i)+1;
                    intersect(i,3) = X_temp;
                    intersect(i,4) = Y_temp;
                end
            end
            
            X_temp = ceil((longside - b(i))/a(i));
            Y_temp = longside;
            if X_temp>0 && X_temp<=shortside
                if l(i) == 0
                    l(i) = 1;
                    intersect(i,1) = X_temp;
                    intersect(i,2) = Y_temp;
                else
                    l(i) = l(i)+1;
                    intersect(i,3) = X_temp;
                    intersect(i,4) = Y_temp;
                end
            end
            
            X_temp = shortside;
            Y_temp = a(i)*shortside +b(i);
            if Y_temp>0 && Y_temp<=longside
                if l(i) == 0
                    l(i) = 1;
                    intersect(i,1) = X_temp;
                    intersect(i,2) = Y_temp;
                else
                    l(i) = l(i)+1;
                    intersect(i,3) = X_temp;
                    intersect(i,4) = Y_temp;
                end
            end
        end
        
        %% Draw lines
        image = zeros(longside,shortside);
        image_MT_gray = image;
        for i = 1:MTnumber
            image_MT = insertShape(image,'line',intersect(i,:), 'linewidth', 3, 'Color', [I*bundled(i) I*bundled(i) I*bundled(i)]);
            image_MT_gray = image_MT_gray + image_MT(:,:,1);
        end
        
        image_MT_gray = image_MT_gray + 10;
        image_MT_gray(image_MT_gray>255) = 255;
        if G>0
            image_MT_gray = imgaussfilt(image_MT_gray,G);
        end
        
        
        
        %% MTSD
        %Apply Sobel Filter over a MTs image to test it
        clear H_full V_full H V M D x y mxd_thr mxd_corrected mxd_indexed
        object_double = image_MT_gray;
        H_full = conv2(object_double,Gx);
        V_full = conv2(object_double,Gy);
        H = H_full(5:longside,5:shortside);
        V = V_full(5:longside,5:shortside);
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
        
        m_added_norm(:,k+1) = m_added/sum(m_added);
        
        %%Length of microtubules
        for i = 1:MTnumber
            Length(i) = sqrt((intersect(i,1)-intersect(i,3))*(intersect(i,1)-intersect(i,3))...
                + (intersect(i,2)-intersect(i,4))*(intersect(i,2)-intersect(i,4)));
        end
        
        for i=1:(length(binrange)-1)
            Value(i,k+1)=0;
            for l=1:MTnumber
                if angles(l)>=binrange(i) && angles(l)<binrange(i+1)
                    Value(i,k+1) = Value(i,k+1) + Length(l)*bundled(l);
                end
            end
        end
        Value(:,k+1) = Value(:,k+1) / sum(Value(:,k+1));
    end
    
    cd(currdir);
    %% Von Mises
    vonmises_fit_dist_sum;
    %% Summary
    summary = zeros(length(SD'),5);
    for counter2 = 1:length(SD')
        summary(counter2,1) = counter2;
    end
    % MTSD
    summary(:,2) = SD';
    % MT direction
    mu(mu<0) = 180+mu(mu<0);
    summary(:,3) = mu';
    m_added_norm = Value;
    vonmises_fit_dist_sum;
    % MTSD theoretical
    summary(:,4) = SD';
    % MT theoretical
    summary(:,5) = mu'+90;
    
    headers2 = {'Cell', 'MTSD', 'MT direction', 'MTSD theor', 'MT direction theor'};
    summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
        num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_',num2str(method), '_' num2str(G),'.csv'];
    cd(result_dir);
    csvwrite_with_headers(summary_filename,summary,headers2);
    
    % MTSD
    Averages(1,1) = mean(summary(:,2));
    Averages(1,2) = sqrt(var(summary(:,2))/length(summary(:,2)));
    % MT direction
    Averages(1,3) = mean(summary(:,3));
    Averages(1,4) = sqrt(var(summary(:,3))/length(summary(:,3)));
    % MTSD theor
    Averages(1,5) = mean(summary(:,4));
    Averages(1,6) = sqrt(var(summary(:,4))/length(summary(:,4)));
    % MT direction theor
    Averages(1,7) = mean(summary(:,5));
    Averages(1,8) = sqrt(var(summary(:,5))/length(summary(:,5)));
    headers = {'MTSD','sem', 'MT direction','sem', ...
        'MTSD theor','sem', 'MT direction theor','sem'};
    summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
        num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_', num2str(method),'_' num2str(G),'_summary.csv'];
   csvwrite_with_headers(summary_filename,Averages,headers);
    cd(currdir);
end


cd(currdir);
clear variables
close all
clc