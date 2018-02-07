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
%% Setting parameters
m_added_norm = zeros(45,cells+1);
m_added_norm(:,1) = bincenter;
Value = zeros(45,cells+1);
Value(:,1) = bincenter;
mts_density = zeros(1, cells);
mts_area = zeros(1, cells);
Uniformity = zeros(1, cells);
Spars = zeros(1, cells);
mts_bundling = zeros(1, cells);
kurt = zeros(1, cells);
skew = zeros(1, cells);
threshold = zeros(1, cells);
Length = zeros(1,MTnumber);



result_dir = '/Users/nataliabulgakova/Desktop/test';
cd('/Users/nataliabulgakova/Desktop/test');
% image_dir_name = ['MTs_', num2str(MTnumber), '_SD' num2str(distribution),'_int',...
%     num2str(I),'_bund',num2str(bundling),'_Ecc',num2str(Eccentricity)];
% if exist([result_dir,'/', image_dir_name],'dir') == 0
%     mkdir(result_dir,image_dir_name);
% end
% image_dir = [result_dir,'/', image_dir_name];
% cd(image_dir);