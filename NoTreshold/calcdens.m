    %% Density data from signal and background
    selected_signal = ones(longside,shortside);
    to_analyse_all = regionprops(selected_signal, image_MT_gray,'PixelValues');   
    [px, py] = gradient(image_MT_gray);
    % Uniformity
    Uniformity(k) = 100 * (1 - sqrt(var(to_analyse_all.PixelValues))/mean(to_analyse_all.PixelValues));
    %Sparseness
    Spars(k) = calcSparseness(to_analyse_all.PixelValues/mean(to_analyse_all.PixelValues),1);
    
    % Kurtosis and skewness
    
    kurt(k) = kurtosis(to_analyse_all.PixelValues);
    skew(k) = skewness(to_analyse_all.PixelValues);
    
    %Sdr - the developed interfacial area ratio
   
    Sdr(k) = (sqrt(1+sum(px(:).*px(:)+(py(:).*py(:))))-1)/length(to_analyse_all.PixelValues);
    %Sdq - the root mean square gradient
    Sdq(k) = sqrt(sum(px(:).*px(:)+(py(:).*py(:)))/length(to_analyse_all.PixelValues));
    
    % Second option
    SdrM(k) = Sdr(k)/(max(to_analyse_all.PixelValues)-min(to_analyse_all.PixelValues));
    %Sdq - the root mean square gradient
    SdqM(k) = Sdq(k)/(max(to_analyse_all.PixelValues)-min(to_analyse_all.PixelValues));
    %Sal - the autocorrelation length
%     C2 = abs(image_MT_gray(2:end,2:end) - image_MT_gray(1:end-1,1:end-1))/mean(to_analyse_all.PixelValues);
%     C = xcorr2(C2);
    close all