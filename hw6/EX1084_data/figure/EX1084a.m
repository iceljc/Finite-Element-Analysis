clc;clear all;close all;
%%
if 1
f1 = load('EX1084_data\U_RE0L.txt');
f2 = load('EX1084_data\U_RE500L.txt');
f3 = load('EX1084_data\U_RE1000L.txt');

x = f1(:,2);
y = f1(:,3);
vx1 = f1(:,4);
vx2 = f2(:,4);
vx3 = f3(:,4);

y_plot1 = y(find(x==0.5));
vx1_plot = vx1(find(x==0.5));
vx2_plot = vx2(find(x==0.5));
vx3_plot = vx3(find(x==0.5));

f4 = load('EX1084_data\U_RE0Q.txt');
f5 = load('EX1084_data\U_RE500Q.txt');
f6 = load('EX1084_data\U_RE1000Q.txt');

x = f4(:,2);
y = f4(:,3);
vx4 = f4(:,4);
vx5 = f5(:,4);
vx6 = f6(:,4);

y_plot2 = y(find(x==0.5));
vx4_plot = vx4(find(x==0.5));
vx5_plot = vx5(find(x==0.5));
vx6_plot = vx6(find(x==0.5));


figure
plot(vx1_plot, y_plot1, '-ok', 'LineWidth', 2, 'MarkerSize', 10); grid on;hold on
plot(vx2_plot, y_plot1, '-sk', 'LineWidth', 2, 'MarkerSize', 10); 
plot(vx3_plot, y_plot1, '-dk', 'LineWidth', 2, 'MarkerSize', 10); 

plot(vx4_plot, y_plot2, '--ok', 'LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0, 0, 0]); 
plot(vx5_plot, y_plot2, '--sk', 'LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0, 0, 0]); 
plot(vx6_plot, y_plot2, '--dk', 'LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0, 0, 0]); hold off;

xlabel('Horizontal velocity, u(0.5,y)', 'FontSize', 20);
ylabel('Distance, y', 'FontSize', 20);
legend('Re = 0, 16x20L', 'Re = 500, 16x20L', 'Re = 1000, 16x20L', 'Re = 0, 8x10Q9', ...
     'Re = 500, 8x10Q9',  'Re = 1000, 8x10Q9', 'Location', 'SouthEast');
axis([-0.4, 1, 0, 1])
% set(gcf, 'position', [10,10,500,500]);
set(gca, 'xtick', -0.4:0.2:1, 'ytick', 0:0.1:1)
set(gca, 'FontSize', 20)
end

%%
if 0
f1 = load('EX1084_data\P_RE0L.txt');
f2 = load('EX1084_data\P_RE500L.txt');
f3 = load('EX1084_data\P_RE1000L.txt');

x = f1(:,1);
y = f1(:,2);
y_top1 = max(y);
pres1 = f1(:, end);
pres2 = f2(:, end);
pres3 = f3(:, end);

x_plot1 = x(find(y==y_top1));
pres1_plot = pres1(find(y==y_top1));
pres2_plot = pres2(find(y==y_top1));
pres3_plot = pres3(find(y==y_top1));

f4 = load('EX1084_data\P_RE0Q.txt');
f5 = load('EX1084_data\P_RE500Q.txt');
f6 = load('EX1084_data\P_RE1000Q.txt');

x = f4(:,1);
y = f4(:,2);
y_top2 = max(y);

pres4 = f4(:, end);
pres5 = f5(:, end);
pres6 = f6(:, end);

x_plot2 = x(find(y==y_top2));
pres4_plot = pres4(find(y==y_top2));
pres5_plot = pres5(find(y==y_top2));
pres6_plot = pres6(find(y==y_top2));

figure
% plot(x_plot, pres1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); 
plot(x_plot1, pres2_plot, '-sk', 'LineWidth', 2, 'MarkerSize', 10); grid on;hold on
plot(x_plot1, pres3_plot, '-dk', 'LineWidth', 2, 'MarkerSize', 10); 

% plot(pres4_plot, x_plot, '--ok', 'LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0, 0, 0]); 
plot(x_plot2, pres5_plot, '--sk', 'LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0, 0, 0]); 
plot(x_plot2, pres6_plot, '--dk', 'LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0, 0, 0]); hold off;

xlabel('Horizontal distance, x', 'FontSize', 20);
ylabel('Pressure, P', 'FontSize', 20);
legend('Re = 500, 16x20L', 'Re = 1000, 16x20L','Re = 500, 8x10Q9',  'Re = 1000, 8x10Q9', 'Location', 'NorthWest');
axis([0, 1, -40, 360])
% set(gcf, 'position', [10,10,500,500]);
set(gca, 'xtick', 0:0.1:1, 'ytick', -40:40:360)
set(gca, 'FontSize', 20)
end












