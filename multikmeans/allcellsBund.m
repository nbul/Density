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

Averages(counterB,1) = bundling;
% Uniformity
Averages(counterB,2) = mean(summary(:,2));
Averages(counterB,3) = sqrt(var(summary(:,2))/length(summary(:,2)));
% Sparseness
Averages(counterB,4) = mean(summary(:,3));
Averages(counterB,5) = sqrt(var(summary(:,3))/length(summary(:,3)));
% Skewness
Averages(counterB,6) = mean(summary(:, 4));
Averages(counterB,7) = sqrt(var(summary(:, 4))/length(summary(:, 4)));
% Kurtosis
Averages(counterB,8) = mean(summary(:, 5));
Averages(counterB,9) = sqrt(var(summary(:, 5))/length(summary(:, 5)));
% MTSD
Averages(counterB,10) = mean(summary(:, 6));
Averages(counterB,11) = sqrt(var(summary(:, 6))/length(summary(:, 6)));
% MT direction
Averages(counterB,12) = mean(summary(:, 7));
Averages(counterB,13) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
% MTSD theoretical
Averages(counterB,14) = mean(summary(:, 8));
Averages(counterB,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));
% MT direction theoretical
Averages(counterB,16) = mean(summary(:, 9));
Averages(counterB,17) = sqrt(var(summary(:, 9))/length(summary(:, 9)));

% Sdr
Averages(counterB,18) = mean(summary(:,10));
Averages(counterB,19) = sqrt(var(summary(:,10))/length(summary(:,10)));
% Sdq
Averages(counterB,20) = mean(summary(:,11));
Averages(counterB,21) = sqrt(var(summary(:,11))/length(summary(:,11)));
% Threshold
Averages(counterB,22) = mean(summary(:,12));
Averages(counterB,23) = sqrt(var(summary(:,12))/length(summary(:,12)));
% Area
Averages(counterB,24) = mean(summary(:,13));
Averages(counterB,25) = sqrt(var(summary(:,13))/length(summary(:,13)));
% Bundling
Averages(counterB,26) = mean(summary(:,14));
Averages(counterB,27) = sqrt(var(summary(:,14))/length(summary(:,14)));

Averages(counterB,28) = length(summary(:,1));