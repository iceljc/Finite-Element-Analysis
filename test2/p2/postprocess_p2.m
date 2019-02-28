clc;clear all;close all;
%%
f1 = load('LOAD-CENTER-W-CCFF.dat');
f2 = load('LOAD-CENTER-W-FFCC.dat');

f3 = load('LOAD-BOUND-W-CCFF.dat');
f4 = load('LOAD-BOUND-W-FFCC.dat');

load_P = f1(:,1)*10;

center_defln_ccff = f1(:,2);
center_defln_ffcc = f2(:,2);

bound_defln_ccff = f3(:,2);
bound_defln_ffcc = f4(:,2);

figure
plot(load_P, center_defln_ccff, '-k', 'LineWidth', 2);grid on; hold on;
plot(load_P, center_defln_ffcc, '--k', 'LineWidth', 2);
plot(load_P, bound_defln_ccff, 'ok', 'LineWidth', 2, 'MarkerSize', 10);
plot(load_P, bound_defln_ffcc, 'sk', 'LineWidth', 2, 'MarkerSize', 10);hold off;

legend('w(0,0), CCFF','w(0,0), FFCC', 'w(a/2,0), CCFF', 'w(0,a/2), FFCC','Location', 'NorthWest');

xlabel('Load Parameter, P', 'FontSize', 20);
ylabel('Displacement, w', 'FontSize', 20);

set(gca, 'FontSize', 20)

%%
f5 = load('STRESS-CCFF.txt');
f6 = load('STRESS-FFCC.txt');

strs_xx_ccff = f5(:,3);
strs_yy_ccff = f5(:,4)/10;
strs_xy_ccff = f5(:,5);

strs_xx_ffcc = f6(:,3)/10;
strs_yy_ffcc = f6(:,4);
strs_xy_ffcc = f6(:,5);

figure
plot(load_P, strs_xx_ccff,  '-ok', 'LineWidth', 2, 'MarkerSize', 10);grid on; hold on;
plot(load_P, strs_yy_ccff,  '--sk', 'LineWidth', 2, 'MarkerSize', 10);
plot(load_P, strs_xy_ccff,  '-.*k', 'LineWidth', 2, 'MarkerSize', 10);

h = legend('$$\overline{\sigma}_{xx}$', '$$\overline{\sigma}_{yy}*0.1$', '$$\overline{\sigma}_{xy}$','Location', 'NorthWest');
set(h, 'Interpreter','latex', 'FontSize', 30);

xlabel('Load Parameter, P', 'FontSize', 20);
ylabel('Dimensionless stress', 'FontSize', 20);
set(gca, 'FontSize', 20)

figure
plot(load_P, strs_xx_ffcc,  '-ok', 'LineWidth', 2, 'MarkerSize', 10);grid on; hold on;
plot(load_P, strs_yy_ffcc,  '--sk', 'LineWidth', 2, 'MarkerSize', 10);
plot(load_P, strs_xy_ffcc,  '-.*k', 'LineWidth', 2, 'MarkerSize', 10);

h = legend('$$\overline{\sigma}_{xx}*0.1$', '$$\overline{\sigma}_{yy}$', '$$\overline{\sigma}_{xy}$','Location', 'NorthWest');
set(h, 'Interpreter','latex', 'FontSize', 30);

xlabel('Load Parameter, P', 'FontSize', 20);
ylabel('Dimensionless stress', 'FontSize', 20);
set(gca, 'FontSize', 20)









