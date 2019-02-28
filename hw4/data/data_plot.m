clc;clear all;close all

f = load('final_output.txt');

center_defln_ss3_L88 = f(:, 1);
center_defln_ss3_Q44 = f(:, 4);

max_stress_xx_ss3_L88 = f(:, 2);
max_stress_xx_ss3_Q44 = f(:, 5);

max_stress_xy_ss3_L88 = f(:, 3);
max_stress_xy_ss3_Q44 = f(:, 6);



center_defln_cc_L88 = f(:, 7);
center_defln_cc_Q44 = f(:, 10);

max_stress_xx_cc_L88 = f(:, 8);
max_stress_xx_cc_Q44 = f(:, 11);

max_stress_xy_cc_L88 = f(:, 9);
max_stress_xy_cc_Q44 = f(:, 12);

Load_P = 5:5:160;

%%
figure
plot(Load_P, center_defln_ss3_L88, '-k', 'LineWidth', 2.5); hold on; grid on;
plot(Load_P, center_defln_cc_L88, '--k','LineWidth', 2.5); 

plot(Load_P, center_defln_ss3_Q44, 'ok', 'LineWidth', 2.5, 'MarkerSize', 10);
plot(Load_P, center_defln_cc_Q44, 'sk', 'LineWidth', 2.5,'MarkerSize', 10); hold off;

xlabel('Load Parameter, P', 'FontSize', 25);
ylabel('Center deflection, W/h', 'FontSize', 25);
legend('SS-3 (8\times8L4)', 'Clamped bc (8\times8L4)', 'SS-3 (4\times4Q9)', 'Clamped bc (4\times4Q9)', 'Location', 'NorthWest')
set(gca, 'FontSize', 25)

figure
plot(Load_P, max_stress_xx_ss3_L88, '-k', 'LineWidth', 2.5); hold on; grid on;
plot(Load_P, max_stress_xx_cc_L88, '--k','LineWidth', 2.5); 

plot(Load_P, max_stress_xx_ss3_Q44, 'ok', 'LineWidth', 2.5, 'MarkerSize', 10);
plot(Load_P, max_stress_xx_cc_Q44, 'sk', 'LineWidth', 2.5, 'MarkerSize', 10); hold off;

xlabel('Load Parameter, P', 'FontSize', 25);
ylabel('Maximum stress, $\overline {\sigma}_{xx,max}$', 'interpreter', 'latex', 'FontSize', 25);
legend('SS-3 (8\times8L4)', 'Clamped bc (8\times8L4)', 'SS-3 (4\times4Q9)', 'Clamped bc (4\times4Q9)', 'Location', 'NorthWest')
set(gca, 'FontSize', 25)

figure
plot(Load_P, max_stress_xy_ss3_L88, '-k', 'LineWidth', 2.5); hold on; grid on;
plot(Load_P, max_stress_xy_cc_L88, '--k', 'LineWidth', 2.5); 

plot(Load_P, max_stress_xy_ss3_Q44, 'ok', 'LineWidth', 2.5, 'MarkerSize', 10); 
plot(Load_P, max_stress_xy_cc_Q44, 'sk', 'LineWidth', 2.5, 'MarkerSize', 10); hold off;

xlabel('Load Parameter, P', 'FontSize', 25);
ylabel('Maximum stress, $\overline {\sigma}_{xy,max}$', 'interpreter', 'latex', 'FontSize', 25);
legend('SS-3 (8\times8L4)', 'Clamped bc (8\times8L4)', 'SS-3 (4\times4Q9)', 'Clamped bc (4\times4Q9)', 'Location', 'NorthWest')
set(gca, 'FontSize', 25)

%%
% figure
% plot(Load_P, center_defln_ss3_Q44, '-k', 'LineWidth', 3); hold on; grid on;
% plot(Load_P, center_defln_cc_Q44, '--k','LineWidth', 3); hold off;
% xlabel('Load Parameter, P', 'FontSize', 25);
% ylabel('Center deflection, W/h', 'FontSize', 25);
% legend('SS-3 (4\times4Q9)', 'Clamped bc (4\times4Q9)', 'Location', 'NorthWest')
% set(gca, 'FontSize', 25)
% 
% figure
% plot(Load_P, max_stress_xx_ss3_Q44, '-k', 'LineWidth', 3); hold on; grid on;
% plot(Load_P, max_stress_xx_cc_Q44, '--k','LineWidth', 3); hold off;
% xlabel('Load Parameter, P', 'FontSize', 25);
% ylabel('Maximum stress, $\overline {\sigma}_{xx,max}$', 'interpreter', 'latex', 'FontSize', 25);
% legend('SS-3 (4\times4Q9)', 'Clamped bc (4\times4Q9)', 'Location', 'NorthWest')
% set(gca, 'FontSize', 25)
% 
% figure
% plot(Load_P, max_stress_xy_ss3_Q44, '-k', 'LineWidth', 3); hold on; grid on;
% plot(Load_P, max_stress_xy_cc_Q44, '--k','LineWidth', 3); hold off;
% xlabel('Load Parameter, P', 'FontSize', 25);
% ylabel('Maximum stress, $\overline {\sigma}_{xy,max}$', 'interpreter', 'latex', 'FontSize', 25);
% legend('SS-3 (4\times4Q9)', 'Clamped bc (4\times4Q9)', 'Location', 'NorthWest')
% set(gca, 'FontSize', 25)








