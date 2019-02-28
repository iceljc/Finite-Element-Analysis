clc;clear all;close all;

f1 = load('P1DILIN.txt');
f2 = load('P1DINLIN.txt');
f3 = load('P1NILIN.txt');
f4 = load('P1NINLIN.txt');

x = f1(1:9, 2);

T_di_lin = f1(1:9, end);
T_di_nlin = f2(1:9, end);
T_ni_lin = f3(1:9, end);
T_ni_nlin = f4(1:9, end);

figure
plot(x, T_di_lin, '-k', 'LineWidth', 2); grid on; hold on;
plot(x, T_di_nlin, '--k', 'LineWidth', 2); 
plot(x, T_ni_lin, 'ok', 'LineWidth', 2, 'MarkerSize', 10); 
plot(x, T_ni_nlin, 'sk', 'LineWidth', 2, 'MarkerSize', 10); hold off;

xlabel('X', 'FontSize', 20);
ylabel('T(x,0)', 'FontSize', 20);
legend('Linear solution with direct iteration', 'Nonlinear solution with direct iteration',...
    'Linear solution with Newton iteration', 'Nonlinear solution with Newton iteration');
set(gca, 'FontSize', 20)











