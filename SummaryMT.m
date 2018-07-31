%% Summary

SDexp = SD;
mu(mu<0) = 180+mu(mu<0);
muexp = mu;
m_added_norm = Value;
vonmises;
mu = mu + 90;
mu(mu<0) = 180+mu(mu<0);

summaryMT = [(1:1:cells)',intensity',intensity_mts', mts_area', mts_density'...
    mts_bundling', Uniformity', UNAAD', Spars', skew', kurt', space',...
    Ent1',Ent2', Sdr', Sdq', SDexp', muexp', SD', mu', mkthr', levelsmk',...
    otsuthr', thr2kmeans', thredges', thrmulti', levelsmulti'];

%% Averaged values
Averages_MT(counterMT,1) = MTnumber;
Averages_MT(counterMT,2:2:size(summaryMT(:,2:end),2)*2) = mean(summaryMT(:,2:end),1);
Averages_MT(counterMT,3:2:(size(summaryMT(:,2:end),2)*2+1))...
    = sqrt(var(summaryMT(:,2:end),1)/size(summaryMT(:,2:end),1));
Averages_MT(counterMT,size(summaryMT(:,2:end),2)*2+2) = cells;


