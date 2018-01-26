
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



result_dir = '/data/md1nbu/Density/Validation';
cd('/data/md1nbu/Density/Validation');
