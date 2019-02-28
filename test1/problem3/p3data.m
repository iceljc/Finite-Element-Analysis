clc;clear all

f = load('p3data.txt');

P = f(:,1);
w0 = f(:,2);


plot(P, w0, '-ok', 'LineWidth', 3);grid on

xlabel('Load parameter, P', 'FontSize',25);
ylabel('Center deflection, W(0)', 'FontSize',25);
set(gca, 'FontSize', 25)
