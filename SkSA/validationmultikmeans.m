clear variables
clc
currdir = pwd;
addpath(pwd);
warning('off','stats:kmeans:FailedToConvergeRep');
cells = 10;
I = 25;
setup;
Eccentricity = 0.8;
setcell;

averages = zeros(10,25);
averages(2:end,1) = (1:1:9) * 5 - 5;
averages(1,2:end) = (1:1:24) * 10 + 10;
averagesSq = averages;
averagesSqsem = averages;
averagesSA = averages;
averagesSAsem = averages;

for Number1 = 1:1:24
    MTnumber = Number1*10 + 10;
    for Number2 = 1:1:9   
        bundling = Number2 * 5 - 5;
        for k=1:cells
            MT;
            calcdens;
        end
        cd(currdir);
        averagesSq(Number2+1,Number1+1) = mean(skew);
        averagesSqsem(Number2+1,Number1+1) = sqrt(var(skew)/length(skew));
        averagesSA(Number2+1,Number1+1) = mean(mts_area);
        averagesSAsem(Number2+1,Number1+1) = sqrt(var(mts_area)/length(mts_area));
    end
end


cd(result_dir);
csvwrite('Skewness.csv',AveragesSq);
csvwrite('Skewness_sem.csv',AveragesSqsem);
csvwrite('SignalArea.csv',AveragesSA);
csvwrite('SignalArea_sem.csv',AveragesSAsem);
cd(currdir);


clear variables
close all
clc