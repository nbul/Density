shortside = 200;
long_side = zeros(1, (4*shortside));
for i=(2*shortside):(10*shortside)
    selected_signal = ones(i,shortside);
    to_analyse_o = regionprops(selected_signal, selected_signal,'Eccentricity');
    long_side(i-2*shortside+1) = abs(to_analyse_o.Eccentricity-Eccentricity);
end
[Ecc, Ind] = min(long_side);
longside = Ind + 2*shortside-1;
