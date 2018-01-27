% function fig_to_lcn(A,xn);
% A = the grey chanel of the figure. 
%      to get the grey chanel out of the image 'name_figure' do
%       %
%       %[x1,map] = imread(name_figure);  % read an image
%       %x1 = double(x1);
%       %A = x1(:,:,1);     % take only the grey channel
%
function V = fig_to_LCN(A,xn)

    [x,y] = meshgrid([-xn:xn]);
    s_average    = exp(-(x/xn).^2  - (y/xn).^2);
    s_average = s_average/sum(sum((s_average)));
 
    D = A - conv2(A,s_average,'same');
    S = sqrt(conv2(D.^2, s_average,'same'));
    c = mean(mean(S));
    V =  D./max(c,S);
    V = V- min(min(V));

end
