%% Summary
summary = zeros(length(mts_area),8);
for counter2 = 1:length(mts_area)
    summary(counter2,1) = counter2;
end
% Signal area
summary(:,2) = mts_area';
% Density
summary(:,3) = mts_density';
% Bundling
summary(:,4) = mts_bundling';
% Uniformity
summary(:,5) = Uniformity';
% Sparseness
summary(:,6) = Spars';
% Skewness
summary(:,7) = skew';
% Kurtosis
summary(:,8) = kurt';
% MTSD
summary(:,9) = SD';
% MT direction
mu(mu<0) = 180+mu(mu<0);
summary(:,10) = mu';
m_added_norm = Value;
vonmises_fit_dist_sum;
% MTSD theoretical
summary(:,11) = SD';
% MT theoretical
summary(:,12) = mu'+90;
summary(:,13) = threshold';

headers2 = {'Cell', 'Signal area', 'Density','Bundling', 'Uniformity', ...
    'Sparseness', 'Skewness', 'Kurtosis', 'MTSD', 'MT direction', 'MTSD theor', 'MT direction theor', 'threshold'};
summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'.csv'];
cd(result_dir);
csvwrite_with_headers(summary_filename,summary,headers2);
cd(currdir);
% Signal area
Averages(1,1) = mean(summary(:,2));
Averages(1,2) = sqrt(var(summary(:,2))/length(summary(:,2)));
% Density
Averages(1,3) = mean(summary(:,3));
Averages(1,4) = sqrt(var(summary(:,3))/length(summary(:,3)));
% Bundling
Averages(1,5) = mean(summary(:,4));
Averages(1,6) = sqrt(var(summary(:,4))/length(summary(:,4)));
% Uniformity
Averages(1,7) = mean(summary(:,5));
Averages(1,8) = sqrt(var(summary(:,5))/length(summary(:,5)));
% Sparseness
Averages(1,9) = mean(summary(:,6));
Averages(1,10) = sqrt(var(summary(:,6))/length(summary(:,6)));
% Skewness
Averages(1,11) = mean(summary(:, 7));
Averages(1,12) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
% Kurtosis
Averages(1,13) = mean(summary(:, 8));
Averages(1,14) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
% MTSD
Averages(1,15) = mean(summary(:, 9));
Averages(1,16) = sqrt(var(summary(:, 9))/length(summary(:, 9)));
% MT direction
Averages(1,17) = mean(summary(:, 10));
Averages(1,18) = sqrt(var(summary(:, 10))/length(summary(:, 10)));
% MTSD theoretical
Averages(1,19) = mean(summary(:, 11));
Averages(1,20) = sqrt(var(summary(:, 11))/length(summary(:, 11)));
% MT direction theoretical
Averages(1,21) = mean(summary(:, 12));
Averages(1,22) = sqrt(var(summary(:, 12))/length(summary(:, 12)));
Averages(1,23) = length(summary(:,1));
headers = {'Signal area', 'sem','Density','sem','Bundling','sem', 'Uniformity','sem', ...
    'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
    'MTSD','sem', 'MT direction','sem', ...
    'MTSD theor','sem', 'MT direction theor','sem', 'Cell number'};
summary_filename = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
    num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity),'_summary.csv'];

cd(result_dir);
csvwrite_with_headers(summary_filename,Averages,headers);
cd(currdir);