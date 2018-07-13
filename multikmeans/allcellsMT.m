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
% MTSD
summary(:,5) = SD';
% MT direction
mu(mu<0) = 180+mu(mu<0);
summary(:,6) = mu';
m_added_norm = Value;
vonmises_fit_dist_sum;
% MTSD theoretical
summary(:,7) = SD';
% MT theoretical
summary(:,8) = mu'+90;
% Sdr
summary(:,9) = Sdr';
% Sdq
summary(:,10) = Sdq';
% Threshold 
summary(:,11) = threshold';
% Area
summary(:,12) = mts_area';
% Bundling
summary(:,13) = mts_bundling';
summary(:,14) = Nculsters';
summary(:,15) = space';
summary(:,16) = Skew2';
summary(:,17) = Ent1';
summary(:,18) = Ent2';

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
% MTSD
Averages(counterMT,8) = mean(summary(:, 5));
Averages(counterMT,9) = sqrt(var(summary(:, 5))/length(summary(:, 5)));
% MT direction
Averages(counterMT,10) = mean(summary(:, 6));
Averages(counterMT,11) = sqrt(var(summary(:, 6))/length(summary(:, 6)));
% MTSD theoretical
Averages(counterMT,12) = mean(summary(:, 7));
Averages(counterMT,13) = sqrt(var(summary(:, 7))/length(summary(:, 7)));
% MT direction theoretical
Averages(counterMT,14) = mean(summary(:, 8));
Averages(counterMT,15) = sqrt(var(summary(:, 8))/length(summary(:, 8)));

% Sdr
Averages(counterMT,16) = mean(summary(:,9));
Averages(counterMT,17) = sqrt(var(summary(:,9))/length(summary(:,9)));
% Sdq
Averages(counterMT,18) = mean(summary(:,10));
Averages(counterMT,19) = sqrt(var(summary(:,10))/length(summary(:,10)));
% Threshold
Averages(counterMT,20) = mean(summary(:,11));
Averages(counterMT,21) = sqrt(var(summary(:,11))/length(summary(:,11)));
% Area
Averages(counterMT,22) = mean(summary(:,12));
Averages(counterMT,23) = sqrt(var(summary(:,12))/length(summary(:,12)));
% Bundling
Averages(counterMT,24) = mean(summary(:,13));
Averages(counterMT,25) = sqrt(var(summary(:,13))/length(summary(:,13)));
%N clusters
Averages(counterMT,26) = mean(summary(:,14));
Averages(counterMT,27) = sqrt(var(summary(:,14))/length(summary(:,14)));
%Space
Averages(counterMT,28) = mean(summary(:,15));
Averages(counterMT,29) = sqrt(var(summary(:,15))/length(summary(:,15)));
%Skewness 2
Averages(counterMT,30) = mean(summary(:,16));
Averages(counterMT,31) = sqrt(var(summary(:,16))/length(summary(:,16)));
% Entropy all
Averages(counterMT,32) = mean(summary(:,17));
Averages(counterMT,33) = sqrt(var(summary(:,17))/length(summary(:,17)));
% Entropy above threshold
Averages(counterMT,34) = mean(summary(:,18));
Averages(counterMT,35) = sqrt(var(summary(:,18))/length(summary(:,18)));
Averages(counterMT,36) = length(summary(:,1));