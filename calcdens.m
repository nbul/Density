    %% Generate Cell Masks.
    signal_original = image_MT_gray .* im_bin_c;
    background_original = image_MT_gray .* (ones(longside,shortside) - im_bin_c);
    
    %% Density data from signal and background
    selected_signal = ones(longside,shortside);
    to_analyse_o = regionprops(selected_signal, signal_original,'PixelValues');
    to_analyse_c = regionprops(selected_signal, im_bin_c,'PixelValues');
    to_analyse_back_o = regionprops(selected_signal, background_original,'PixelValues');
    to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');
    
    % Relative to background and signal area
    sum_pixvalues_o = sum(to_analyse_o.PixelValues(:,1));
    sum_pixvalues_back_o = sum(to_analyse_back_o.PixelValues(:,1));
    num_pixvalues_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) ~= 0,1));
    num_pixvalues_back_c = length(to_analyse_c.PixelValues(to_analyse_c.PixelValues(:, 1) == 0,1));
    
    mts_density(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        (sum_pixvalues_back_o / num_pixvalues_back_c)) * (num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c));
    
    % Signal Area
    mts_area(k) = num_pixvalues_c / (num_pixvalues_c + num_pixvalues_back_c);
    
    if max(to_analyse_c.PixelValues)~= 0
        mts_bundling(k) = (mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))-...
            min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))) / ...
            min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
        %         mts_bundling(k) = (((sum_pixvalues_o / num_pixvalues_c) - (sum_pixvalues_back_o / num_pixvalues_back_c)) / ...
        %             (sum_pixvalues_back_o / num_pixvalues_back_c));
        %         mts_bundling(k) = mean(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1))...
        %             /min(to_analyse_o.PixelValues(to_analyse_c.PixelValues(:,1)~= 0,1));
    else
        mts_bundling(k) = 0;
    end
    
    % Uniformity
    Uniformity(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))./...
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
%     image = figure;
%     imshow(image_MT_gray, [0 255]);
%     image_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',num2str(I),...
%         '_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_', num2str(k), '.tif'];
%     print(image, '-dtiff', '-r150', image_filename);
    close all