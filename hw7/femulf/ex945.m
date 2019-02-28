clc;clear all;

%%
fid1 = load('DISP_5Q8_3X3_WI.dat');
fid2 = load('DISP_5Q8_3X3_WOI.dat');

f0 = fid1(:,end)*(-50);

disp_5q8_3x3_wi = [f0, fid1(:, 2), -fid1(:,3:end-1)];
disp_5q8_3x3_woi = [f0, fid2(:, 2), -fid2(:,3:end-1)];

%%
fid3 = load('DISP_5Q8_2X2_WI.dat');
fid4 = load('DISP_5Q8_2X2_WOI.dat');

disp_5q8_2x2_wi = [f0, fid3(:, 2), -fid3(:,3:end-1)];
disp_5q8_2x2_woi = [f0, fid4(:, 2), -fid4(:,3:end-1)];

%%
fid5 = load('STRESS_5Q8_3X3.dat');
fid6 = load('STRESS_5Q8_2X2.dat');

stress_5q8_3x3 = [f0, fid5(:,1:2), fid5(:,3:4)*1e-5, -fid5(:,5)*1e-5];
stress_5q8_2x2 = [f0, fid6(:,1:2), fid6(:,3:4)*1e-5, -fid6(:,5)*1e-5];

fid7 = load('PKSTRESS_5Q8_3X3.dat');
fid8 = load('PKSTRESS_5Q8_2X2.dat');

pkstress_5q8_3x3 = [f0, fid7(:,1:2), fid7(:,3:4)*1e-5, -fid7(:,end)*1e-5];
pkstress_5q8_2x2 = [f0, fid8(:,1:2), fid8(:,3:4)*1e-5, -fid8(:,end)*1e-5];

%%
fid9 = load('DISP_5Q9_3X3.dat');
f02 = fid9(:,end)*(-50);
disp_5q9_3x3_wi = [f02, fid9(:, 2), -fid9(:, 3:end-1)];

fid10 = load('STRESS_5Q9_3X3.dat');
fid11 = load('PKSTRESS_5Q9_3X3.dat');

stress_5q9_3x3 = [f02, fid10(:, 1:2), fid10(:, 3)*1e-5, -fid10(:, end)*1e-5];
pkstress_5q9_3x3 = [f02, fid11(:, 1:2), fid11(:, 3)*1e-5, -fid11(:, end)*1e-5];

%%
figure
plot(f0, disp_5q8_3x3_wi(:,end-1), '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on;
plot(f0, disp_5q8_3x3_wi(:,end), '-ok', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0,0,0]);
plot(f0, disp_5q8_3x3_woi(:,end-1), '-sk', 'LineWidth', 2, 'MarkerSize', 10);
plot(f0, disp_5q8_3x3_woi(:,end), '-sk', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0,0,0]);
plot(f0, disp_5q8_2x2_wi(:,end-1), '--dk', 'LineWidth', 2, 'MarkerSize', 10);
plot(f0, disp_5q8_2x2_wi(:,end), '--dk', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0,0,0]); hold off
axis([0, 1600, 0, 12])
xlabel('Load, f_{0}', 'FontSize', 20);
ylabel('Displacements, -u, -v', 'FontSize', 20);
legend('-u, with iteration, NGPF=3, NGPR=2', '-v, with iteration, NGPF=3, NGPR=2',...
    '-u, without iteration, NGPF=3, NGPR=2', '-v, without iteration, NGPF=3, NGPR=2',...
    '-u, with iteration, NGPF=2, NGPR=2', '-v, with iteration, NGPF=2, NGPR=2', 'Location', 'NorthWest');
text(1000, 9, 'Vertical displacement, -v', 'FontSize', 20);
text(800, 2, 'Horizontal displacement, -u', 'FontSize', 20)
set(gca, 'FontSize', 20)

figure
plot(f0, stress_5q8_3x3(:, 4), '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on;
plot(f0, stress_5q8_3x3(:, 5), '-sk', 'LineWidth', 2, 'MarkerSize', 10);
plot(f0, stress_5q8_3x3(:, 6), '-dk', 'LineWidth', 2, 'MarkerSize', 10);

plot(f0, stress_5q8_2x2(:, 4), '--ok', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0,0,0]);
plot(f0, stress_5q8_2x2(:, 5), '--sk', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0,0,0]);
plot(f0, stress_5q8_2x2(:, 6), '--dk', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0,0,0]);

plot(f0, pkstress_5q8_3x3(:, 4), '-*k', 'LineWidth', 2, 'MarkerSize', 10); hold off

legend('\sigma_{xx}, NGPF=3, NGPR=2', '\sigma_{yy}, NGPF=3, NGPR=2', '\sigma_{xy}, NGPF=3, NGPR=2',...
    '\sigma_{xx}, NGPF=2, NGPR=2', '\sigma_{yy}, NGPF=2, NGPR=2', '\sigma_{xy}, NGPF=2, NGPR=2', 'S_{xx}, NGPF=3, NGPR=2',...
    'Location', 'NorthWest');

axis([0, 1600, 0, 25])
xlabel('Load, f_{0}', 'FontSize', 20);
ylabel('Stresses', 'FontSize', 20);
set(gca, 'FontSize', 20)























