clc; clear all; close all;

fid1 = load('DISP_FINALP6.dat');
fid2 = load('STRESS_FINALP6.dat');
fid3 = load('PKSTRESS_FINALP6.dat');

load_step = abs(fid1(:, end));

disp_v = -fid1(:, end-1);  disp_y = -fid1(:, 3);
cauchy_sxx = fid2(:, 3) * 1e-5;
pk_sxx = fid3(:, 3) * 1e-5;

figure
plot(load_step, disp_v, '-ok', 'LineWidth', 2, 'MarkerSize', 10); grid on
xlabel('Load step', 'FontSize', 20);
ylabel('Displacement, -v', 'FontSize', 20);
set(gca, 'FontSize', 20)

figure
plot(load_step, cauchy_sxx, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on; grid on
plot(load_step, pk_sxx, '--sk', 'LineWidth', 2, 'MarkerSize', 10); hold off
xlabel('Load step', 'FontSize', 20);
ylabel('Stress', 'FontSize', 20);
legend('Cauchy stress, \sigma_{xx}', 'Second Piola-Kirchhoff stress, S_{xx}', 'Location', 'NorthWest')
set(gca, 'FontSize', 20)





