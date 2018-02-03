%% Summary
summary = zeros(cells,7);
for counter2 = 1:cells
    summary(counter2,1) = counter2;
end
% Uniformity
summary(:,2) = Uniformity';
% Sparseness
summary(:,3) = Spars';
% Skewness
summary(:,4) = skew';
% Kurtosis
summary(:,5) = kurt';
% MTSD
summary(:,6) = SD';
% MT direction
mu(mu<0) = 180+mu(mu<0);
summary(:,7) = mu';
m_added_norm = Value;
vonmises_fit_dist_sum;
% MTSD theoretical
summary(:,8) = SD';
% MT theoretical
summary(:,9) = mu'+90;
% Sdr
summary(:,10) = Sdr';
% Sdq
summary(:,11) = Sdq';
% Threshold 
summary(:,12) = threshold';
% Area
summary(:,13) = mts_area';
% Bundling
summary(:,14) = mts_bundling';

%% Averaged values
Averages(counterMT,1) = MTnumber;
% Uniformity
Averages(counterMT,2) = mean(summary(:,2));
Averages(counterMT,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
% Sparseness
Averages(counterMT,4) = mean(summary(:,3));
Averages(counterMT,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
% Skewness
Averages(counterMT,6) = mean(summary(:, 4));
Averages(counterMT,7) = sqrt(var(summary(:, 4))/length(summary(:, 4)));
% Kurtosis
Averages(counterMT,8) = mean(summary(:, 5));
Averages(counterMT,9) = sqrt(var(summary(:, 5))/length(summary(:, 5)));
% MTSD
Averages(counterMT,10) = mean(summary(:, 6));
Averages(counterMT,11) = sqrt(var(summary(:, 6))/length(summary(:, 6)));
% MT direction
Averages(counterMT,12) = mean(summary(:, 7));
Averages(counterMT,13) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
% MTSD theoretical
Averages(counterMT,14) = mean(summary(:, 8));
Averages(counterMT,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
% MT direction theoretical
Averages(counterMT,16) = mean(summary(:, 9));
Averages(counterMT,17) = sqrt(var(summary(:, 9))/length(summary(:, 9)));

% Sdr
Averages(counterMT,18) = mean(summary(:,10));
Averages(counterMT,19) = sqrt(var(summary(:,10))/length(summary(:,10)));
% Sdq
Averages(counterMT,20) = mean(summary(:,11));
Averages(counterMT,21) = sqrt(var(summary(:,11))/length(summary(:,11)));
% Threshold
Averages(counterMT,22) = mean(summary(:,12));
Averages(counterMT,23) = sqrt(var(summary(:,12))/length(summary(:,12)));
% Area
Averages(counterMT,24) = mean(summary(:,13));
Averages(counterMT,25) = sqrt(var(summary(:,13))/length(summary(:,13)));
% Bundling
Averages(counterMT,26) = mean(summary(:,14));
Averages(counterMT,27) = sqrt(var(summary(:,14))/length(summary(:,14)));

Averages(counterMT,28) = length(summary(:,1));