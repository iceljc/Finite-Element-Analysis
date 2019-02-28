clc;clear all;close all;


fid1 = load('temp_linear.dat');
fid2 = load('temp_nonlin.dat');

u1 = fid1(:,4);
v1 = fid1(:,5);
temp1 = fid1(:,end);

u2 = fid2(:,4);
v2 = fid2(:,5);
temp2 = fid2(:,end);

N = 41;

x = linspace(0, 1, N);
y = linspace(0, 1, N);

[X,Y] = meshgrid(x,y);

T1 = reshape(temp1, N, N);
T2 = reshape(temp2, N, N);
U = reshape(u2, N, N);
V = reshape(v2, N, N);

figure
contourf(X,Y,T1',50,'LineColor', 'none'); colorbar;
caxis([0,1]);
xlabel('X', 'FontSize', 20);
ylabel('Y', 'FontSize', 20);
set(gca, 'FontSize', 20)

figure
contourf(X,Y,T2',50,'LineColor', 'none'); colorbar;
caxis([0,1]);
xlabel('X', 'FontSize', 20);
ylabel('Y', 'FontSize', 20);
set(gca, 'FontSize', 20)

figure
quiver(X,Y,U',V','LineWidth',1,'MaxHeadSize',10,'Color','black')
axis([0,1,0,1]);
xlabel('X', 'FontSize', 20);
ylabel('Y', 'FontSize', 20);
set(gca, 'FontSize', 20)


%%
fid3 = load('temp_nonlin_highra.dat');
u3 = fid3(:,4);
v3 = fid3(:,5);
temp3 = fid3(:,end);

%%%%%%%
x1 = fid1(:, 2);
y1 = fid1(:, 3);

x1_plot = x1(find(y1==0.5));
v1_plot = v1(find(y1==0.5));
T1_ploty = temp1(find(y1==0.5));

y1_plot = y1(find(x1==0.5));
u1_plot = u1(find(x1==0.5));
T1_plotx = temp1(find(x1==0.5));

%%%%%%%%%%%
x2 = fid2(:, 2);
y2 = fid2(:, 3);

x2_plot = x2(find(y2==0.5));
v2_plot = v2(find(y2==0.5));
T2_ploty = temp2(find(y2==0.5));

y2_plot = y2(find(x2==0.5));
u2_plot = u2(find(x2==0.5));
T2_plotx = temp2(find(x2==0.5));

%%%%%%%%%%%%%
x3 = fid3(:, 2);
y3 = fid3(:, 3);

x3_plot = x3(find(y3==0.5));
v3_plot = v3(find(y3==0.5));

y3_plot = y3(find(x3==0.5));
u3_plot = u3(find(x3==0.5));

%%%%%%%%%%%%%%
figure
plot(x1_plot, v1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(x2_plot, v2_plot, '-sk', 'LineWidth', 2, 'MarkerSize', 10);
% plot(x3_plot, v3_plot, '-*k'); 
hold off;
xlabel('x^{*}', 'FontSize', 20);
ylabel('v^{*}', 'FontSize', 20);
legend('Ra = 10^{3}', 'Ra = 10^{4}')%, 'Ra = 6x10^{4}');
set(gca, 'FontSize', 20)

figure
plot(u1_plot, y1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(u2_plot, y2_plot, '-sk', 'LineWidth', 2, 'MarkerSize', 10)
% plot(u3_plot, y3_plot, '-*k'); 
hold off
xlabel('u^{*}', 'FontSize', 20);
ylabel('y^{*}', 'FontSize', 20);
legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Location', 'NorthWest')%, 'Ra = 6x10^{4}');
set(gca, 'FontSize', 20)

figure
plot(x1_plot, T1_ploty, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(x2_plot, T2_ploty, '-sk', 'LineWidth', 2, 'MarkerSize', 10)
hold off
xlabel('x^{*}', 'FontSize', 20);
ylabel('T^{*}', 'FontSize', 20);
axis([0,1,0,1]);
legend('Ra = 10^{3}', 'Ra = 10^{4}')%, 'Ra = 6x10^{4}');
set(gca, 'FontSize', 20)

figure
plot(T1_plotx, y1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(T2_plotx, y2_plot, '-sk', 'LineWidth', 2, 'MarkerSize', 10)
hold off
xlabel('T^{*}', 'FontSize', 20);
ylabel('y^{*}', 'FontSize', 20);
axis([0,1,0,1]);
legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Location', 'NorthWest')%, 'Ra = 6x10^{4}');
set(gca, 'FontSize', 20)









