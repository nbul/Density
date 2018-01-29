    %% Density data from signal and background
    selected_signal = ones(longside,shortside);
    to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');   
    
    % Uniformity
    Uniformity(k) = 100 * (1 - sum(abs(to_analyse_all.PixelValues - mean(to_analyse_all.PixelValues))./...
        mean(to_analyse_all.PixelValues))/length(to_analyse_all.PixelValues));
    %Sparseness
    Spars(k) = calcSparseness(to_analyse_all.PixelValues/mean(to_analyse_all.PixelValues),1);
    
    % Kurtosis and skewness
    
    kurt(k) = kurtosis(to_analyse_all.PixelValues);
    skew(k) = skewness(to_analyse_all.PixelValues);
    close all