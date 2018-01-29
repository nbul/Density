clear variables
clc
currdir = pwd;
addpath(pwd);
warning('off','stats:kmeans:FailedToConvergeRep');
cells = 10;
MTnumber = 25;
I = 25;
setup;

for I = 15:5:40
    for Eccentricity = 0.7:0.05:0.95
        setcell;
        for distribution = 20:5:50
            bundling = 0;
            counterMT = 0;            
            Averages = zeros(1,18);
            for MTnumber = 25:5:250
                counterMT = counterMT +1;
                for k=1:cells
                    MT;
                    calcdens;
                    calcMTSD;
                end
                cd(currdir);
                %% Von Mises
                vonmises_fit_dist_sum;
                allcellsMT;
            end
            headers = {'MT number', 'Uniformity','sem', ...
                'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
                'MTSD','sem', 'MT direction','sem', ...
                'MTSD theor','sem', 'MT direction theor','sem','Cell number'};
            summary_filename = ['MTnumber_','SD', num2str(distribution),'_int',...
                num2str(I),'_Ecc',num2str(Eccentricity),'_summary.csv'];
            
            cd(result_dir);
            csvwrite_with_headers(summary_filename,Averages,headers);
            cd(currdir);
            MTnumber = 100;
            Averages = zeros(1,18);
            counterB = 0;
            for bundling = 0:5:45
                counterB = counterB + 1;
                for k=1:cells
                    MT;
                    % Threshold
                    calcdens;
                    calcMTSD;
                end
                cd(currdir);
                %% Von Mises
                vonmises_fit_dist_sum;
                allcellsBund;
            end
            headers = {'Bundling', 'Uniformity','sem', ...
                'Sparseness','sem', 'Skewness','sem', 'Kurtosis','sem',...
                'MTSD','sem', 'MT direction','sem', ...
                'MTSD theor','sem', 'MT direction theor','sem','Cell number'};
            summary_filename = ['Bundling_','SD', num2str(distribution),'_int',...
                num2str(I),'_Ecc',num2str(Eccentricity),'_summary.csv'];
            cd(result_dir);
            csvwrite_with_headers(summary_filename,Averages,headers);
            cd(currdir);
        end
    end
end
clear variables
close all
clc