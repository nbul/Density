    %% Density data from signal and background
    selected_signal = ones(longside,shortside);
    to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');   
    
    % Uniformity
    Uniformity(k) = 100 * (1 - sqrt(var(to_analyse_all.PixelValues))/mean(to_analyse_all.PixelValues));
    %Sparseness
    Spars(k) = calcSparseness(to_analyse_all.PixelValues/mean(to_analyse_all.PixelValues),1);
    
    % Kurtosis and skewness
    
    kurt(k) = kurtosis(to_analyse_all.PixelValues);
    skew(k) = skewness(to_analyse_all.PixelValues);
    
    %Sdr - the developed interfacial area ratio
   
    Sdr_temp = 0;
    for i = 1:(longside - 1)
        for z = 1:(shortside - 1)
            Sdr_temp = Sdr_temp + (abs(image_MT_gray(i+1,k)-image_MT_gray(i,k))/mean(to_analyse_all.PixelValues))^2;
            Sdr_temp = Sdr_temp + (abs(image_MT_gray(i,k+1)-image_MT_gray(i,k))/mean(to_analyse_all.PixelValues))^2;
        end
    end
    Sdr(k) = (sqrt(Sdr_temp + 1) - 1)/length(to_analyse_all.PixelValues);
    %Sdq - the root mean square gradient
    Sdq(k) = sqrt(Sdr_temp/length(to_analyse_all.PixelValues));
    %Sal - the autocorrelation length
%     C2 = abs(image_MT_gray(2:end,2:end) - image_MT_gray(1:end-1,1:end-1))/mean(to_analyse_all.PixelValues);
%     C = xcorr2(C2);
    close all