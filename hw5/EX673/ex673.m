clc;clear all;close all

ulinCAM = load('ULINCAM.dat');

unlCAM = load('UNLCAM.dat');


time1 = ulinCAM(:,1);
time2 = unlCAM(:,1);

center_disp1 = ulinCAM(:,2);
center_disp2 = unlCAM(:,2);


figure
plot(time1, center_disp1, '-ok', 'LineWidth', 2.5, 'MarkerSize', 10); hold on;grid on
plot(time2, center_disp2, '--sk', 'LineWidth', 2.5, 'MarkerSize', 10); hold off
axis([0,6,-0.5, 1])
xlabel('Time, t', 'FontSize', 20)
ylabel('Displacement, u(0,0,t)', 'FontSize', 20)
legend('Linear case', 'Nonlinear case')
set(gca, 'FontSize', 20)

ulinLAM1 = load('ULINLAM1.dat');
ulinLAM2 = load('ULINLAM2.dat');

