MTnumber = 150;
Eccentricity = 0.92;

cells = 10;
distribution = 30;
I = 25;
bundling = 1;
%% Parameters
bin_size = 4;
binrange = -90 : bin_size : 90;
bincenter=binrange(1:(end-1)) + bin_size/2;
Gx = [-2 -1 0 1 2;-3 -2 0 2 3;-4 -3 0 3 4;-3 -2 0 2 3;-2 -1 0 1 2];
Gy = Gx';
stretch_factor = 0.5;
%% Setting parameters
m_added_norm = zeros(45,cells+1);
m_added_norm(:,1) = bincenter;
Value = zeros(45,cells+1);
Value(:,1) = bincenter;
Length = zeros(1,MTnumber);



result_dir = '/Volumes/DataGurdon/Natalia Bulgakova/MT methods paper/Data simulated';
cd('/Volumes/DataGurdon/Natalia Bulgakova/MT methods paper/Data simulated');

headersMT = {'Number MTs','Cell intensity','sem','MT intensity','sem',...
    'MT area','sem','Density','sem','Bundling','sem','Uniformity','sem',...
    'UNAAD','sem','Sparseness','sem','Skewness','sem','Kurtosis','sem'...
    'Gaps','sem','Entropy','sem','MT Entropy','sem','Sdr','sem','Sdq','sem'...
    'SD', 'sem','DEV','sem','SD theor','sem','DEV theor','sem', ...
    'multi-k-means','sem','Clusters','sem','Otsu','sem','2kmeans','sem'...
    'Edges','sem','Multiotsu','sem', 'Clusters','sem', 'Cells'};
headersB = {'Bundling','Cell intensity','sem','MT intensity','sem',...
    'MT area','sem','Density','sem','Bundling','sem','Uniformity','sem',...
    'UNAAD','sem','Sparseness','sem','Skewness','sem','Kurtosis','sem'...
    'Gaps','sem','Entropy','sem','MT Entropy','sem','Sdr','sem','Sdq','sem'...
    'SD', 'sem','DEV','sem','SD theor','sem','DEV theor','sem', ...
    'multi-k-means','sem','Clusters','sem','Otsu','sem','2kmeans','sem'...
    'Edges','sem','Multiotsu','sem', 'Clusters','sem', 'Cells'};
% image_dir_name = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
%     num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity)];
% if exist([result_dir,'/', image_dir_name],'dir') == 0
%     mkdir(result_dir,image_dir_name);
% end
% image_dir = [result_dir,'/', image_dir_name];
% cd(image_dir);