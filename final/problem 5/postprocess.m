clc; clear all; close all;

%%
fid0 = load('V_RE0.dat');
fid1 = load('V_RE5000.dat');
fid2 = load('V_RE10000.dat');


x = fid1(:, 2);
y = fid1(:, 3);

u0 = fid0(:, 4);
v0 = fid0(:, 5);
u1 = fid1(:, 4);
v1 = fid1(:, 5);
u2 = fid2(:, 4);
v2 = fid2(:, 5);

y_plot = y(find(x == 0.5));
u0_plot = u0(find(x == 0.5));
v0_plot = v0(find(x == 0.5));
u1_plot = u1(find(x == 0.5));
v1_plot = v1(find(x == 0.5));
u2_plot = u2(find(x == 0.5));
v2_plot = v2(find(x == 0.5));

figure
plot(u0_plot, y_plot, '-.*k', 'LineWidth', 2, 'MarkerSize', 10); hold on; grid on
plot(u1_plot, y_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); 
plot(u2_plot, y_plot, '--sk', 'LineWidth', 2, 'MarkerSize', 10); hold off
axis([-0.2, 1, 0, 1])
xlabel('Horizontal velocity, v_{x}', 'FontSize', 20);
ylabel('Distance, y', 'FontSize', 20);
legend('Re = 0', 'Re = 5,000', 'Re = 10,000');
set(gca, 'FontSize', 20);

%%
fid3 = load('P_RE5000.dat');
fid4 = load('P_RE10000.dat');

x = fid3(:, 1);
y = fid3(:, 2);

P1 = fid3(:, end);
P2 = fid4(:, end);

x_plot = x(find(y == 0.15625));
P1_plot = P1(find(y == 0.15625));
P2_plot = P2(find(y == 0.15625));

figure
plot(x_plot, P1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on; grid on
plot(x_plot, P2_plot, '--sk', 'LineWidth', 2, 'MarkerSize', 10); hold off
% axis([-0.2, 1, 0, 1])
xlabel('Distance, x', 'FontSize', 20);
ylabel('Pressure, P', 'FontSize', 20);
legend('Re = 5,000', 'Re = 10,000');
set(gca, 'FontSize', 20);





