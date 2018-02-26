clear variables
clc
currdir = pwd;

Area = csvread('/Volumes/DataGurdon/Natalia Bulgakova/MT methods paper/Fig7/No threshold/SignalArea.csv');
Skewness = csvread('/Volumes/DataGurdon/Natalia Bulgakova/MT methods paper/Fig7/No threshold/Skewness.csv');
A = zeros(216,3);
counter = 0;
ZA = Area(2:end,2:end);
[XA, YA] = meshgrid(Area(1,2:end),Area(2:end,1));
XYA(:,:,1) = XA;
XYA(:,:,2) = YA;
ZS = Skewness(2:end,2:end);
[XS, YS] = meshgrid(Skewness(1,2:end),Skewness(2:end,1));
surf(XA, YA, ZA);
surf(XS, YS, ZS);


fitA = fittype('((a+b*(1-exp(-c*x))))*(d-e*y)',...
    'dependent',{'z'},'independent',{'x','y'},...
    'coefficients',{'a','b','c','d','e'});
fA = fit([XA(:), YA(:)],ZA(:),fitA, 'Lower', [-10 -10 0 -10 -10], 'Upper', [10 20 0.05 20 10], 'Robust', 'LAR');
plot(fA,[XA(:), YA(:)],ZA(:));

fitS = fittype('((a+b*(exp(-c*x))))*(d+e*(1-exp(-f*y)))',...
    'dependent',{'z'},'independent',{'x','y'},...
    'coefficients',{'a','b','c','d','e','f'});
fS = fit([XS(:), YS(:)],ZS(:),fitS, 'Lower', [0 0 0 0 0 0], 'Upper', [1 10 10 10 10 1], 'Robust', 'LAR');
plot(fS,[XS(:), YS(:)],ZS(:));