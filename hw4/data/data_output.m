clc;clear all

center_defln = load('LOAD-CENTER-DEFLE.dat');
stress = load('LOAD-STRESS-MAX.dat');

NL = 32;
FS = stress(:,1);

stress_xx_max = [];
stress_xy_max = [];
stress_xx_max_info = [];
stress_xy_max_info = [];

for i = 1 : NL
   row = find(FS == i); 
   data_mod = stress(row, 2:end);
   [sigmaxx_max, index1] = max(stress(row, 4));
   stress_xx_max_info = [stress_xx_max_info; data_mod(index1, :)];
   stress_xx_max = [stress_xx_max; sigmaxx_max];
   
   [sigmaxy_max, index2] = max(stress(row, end));
   stress_xy_max_info = [stress_xy_max_info; data_mod(index2, :)];
   stress_xy_max = [stress_xy_max; sigmaxy_max];

end
