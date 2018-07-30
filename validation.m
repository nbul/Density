clear variables
clc
currdir = pwd;
addpath(pwd);
warning('off','stats:kmeans:FailedToConvergeRep');
cells = 10;
MTnumber = 25;
I = 25;
setup;


for Eccentricity = 0.7:0.1:0.9
    setcell;
    for I = 15:10:35        
        for distribution = 30:10:40
            bundling = 0;
            counterMT = 0;
            Averages_MT = zeros(1,18);
            for MTnumber = 25:25:250
                counterMT = counterMT +1;
                for k=1:cells
                    MT;
                    threshold;
                    calcdens;
                    calcMTSD;
                end
                cd(currdir);
                %% Von Mises
                vonmises;
                SummaryMT;
            end
           
            summary_filename = ['MTnumber_','SD', num2str(distribution),'_int',...
                num2str(I),'_Ecc',num2str(Eccentricity),'_summary.csv'];
            
            cd(result_dir);
            csvwrite_with_headers(summary_filename,Averages_MT,headersMT);
            cd(currdir);
            MTnumber = 100;
            Averages_B = zeros(1,18);
            counterB = 0;
            for bundling = 0:10:40
                counterB = counterB + 1;
                for k=1:cells
                    MT;
                    threshold;
                    calcdens;
                    calcMTSD;
                end
                cd(currdir);
                %% Von Mises
                vonmises;
                SummaryB;
            end
  
            summary_filename = ['Bundling_','SD', num2str(distribution),'_int',...
                num2str(I),'_Ecc',num2str(Eccentricity),'_summary.csv'];
            cd(result_dir);
            csvwrite_with_headers(summary_filename,Averages_B,headersB);
            cd(currdir);
        end
    end
end
clear variables
close all
clc