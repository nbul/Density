%% Summary

SDexp = SD;
mu(mu<0) = 180+mu(mu<0);
muexp = mu;
m_added_norm = Value;
vonmises;
mu = mu + 90;
mu(mu<0) = 180+mu(mu<0);

summaryB = [(1:1:cells)',intensity',intensity_mts', mts_area', mts_density'...
    mts_bundling', Uniformity', UNAAD', Spars', skew', kurt', space',...
    Ent1',Ent2', Sdr', Sdq', SDexp', muexp', SD', mu', mkthr', levelsmk',...
    otsuthr', thr2kmeans', thredges', thrmulti', levelsmulti'];

%% Averaged values
Averages_B(counterB,1) = bundling;
Averages_B(counterB,2:2:size(summaryB(:,2:end),2)*2) = mean(summaryB(:,2:end),1);
Averages_B(counterB,3:2:(size(summaryB(:,2:end),2)*2+1))...
    = sqrt(var(summaryB(:,2:end),1)/size(summaryB(:,2:end),1));
Averages_B(counterB,size(summaryB(:,2:end),2)*2+2) = cells;


