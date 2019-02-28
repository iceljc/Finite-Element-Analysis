clc;clear all;close all

fid1 = load('CENTER_DEFLN_EBT_N0.dat');
fid2 = load('CENTER_DEFLN_EBT_N1.dat');
fid3 = load('CENTER_DEFLN_EBT_N5.dat');

fid4 = load('CENTER_DEFLN_TBT_N0.dat');
fid5 = load('CENTER_DEFLN_TBT_N1.dat');
fid6 = load('CENTER_DEFLN_TBT_N5.dat');

q = fid1(:,1);

w_center_ebt_n0 = fid1(:, 2);
w_center_ebt_n1 = fid2(:, 2);
w_center_ebt_n5 = fid3(:, 2);

w_center_tbt_n0 = fid4(:, 2);
w_center_tbt_n1 = fid5(:, 2);
w_center_tbt_n5 = fid6(:, 2);

figure
plot(q, w_center_ebt_n0, '-k', 'LineWidth', 2, 'MarkerSize', 15); grid on; hold on
plot(q, w_center_ebt_n1, '--k', 'LineWidth', 2, 'MarkerSize', 15);
plot(q, w_center_ebt_n5, '-.k', 'LineWidth', 2, 'MarkerSize', 15);

plot(q, w_center_tbt_n0, 'ok', 'LineWidth', 2, 'MarkerSize', 15); 
plot(q, w_center_tbt_n1, 'sk', 'LineWidth', 2, 'MarkerSize', 15);
plot(q, w_center_tbt_n5, '*k', 'LineWidth', 2, 'MarkerSize', 15); hold off

xlabel('Load, q0', 'FontSize', 20);
ylabel('Center deflection, w', 'FontSize', 20);
legend('EBT, N=0', 'EBT, N=1', 'EBT, N=5', 'TBT, N=0', 'TBT, N=1', 'TBT, N=5','Location', 'SouthEast');
set(gca, 'FontSize', 20)








