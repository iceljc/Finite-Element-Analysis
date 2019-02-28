clc;clear all;close all;


fid1 = load('temp_ra10e3.dat');
fid2 = load('temp_ra10e4.dat');
fid3 = load('temp_ra10e5.dat');

u1 = fid1(:,4);
v1 = fid1(:,5);
temp1 = fid1(:,end);

u2 = fid2(:,4);
v2 = fid2(:,5);
temp2 = fid2(:,end);

u3 = fid3(:,4);
v3 = fid3(:,5);
temp3 = fid3(:,end);

N = 41;

x = linspace(0, 1, N);
y = linspace(0, 1, N);

[X,Y] = meshgrid(x,y);

T1 = reshape(temp1, N, N);
T2 = reshape(temp2, N, N);
T3 = reshape(temp3, N, N);
U = reshape(u3, N, N);
V = reshape(v3, N, N);

% figure
% contourf(X,Y,T1',50,'LineColor', 'none'); colorbar;
% caxis([0,1]);
% xlabel('X', 'FontSize', 20);
% ylabel('Y', 'FontSize', 20);
% set(gca, 'FontSize', 20)
% 
% figure
% contourf(X,Y,T2',50,'LineColor', 'none'); colorbar;
% caxis([0,1]);
% xlabel('X', 'FontSize', 20);
% ylabel('Y', 'FontSize', 20);
% set(gca, 'FontSize', 20)

% figure
% contourf(X,Y,T3',50,'LineColor', 'none'); colorbar;
% caxis([0,1]);
% xlabel('X', 'FontSize', 20);
% ylabel('Y', 'FontSize', 20);
% set(gca, 'FontSize', 20)
% 
% figure
% quiver(X,Y,U',V','LineWidth',1,'MaxHeadSize',10,'Color','black')
% axis([0,1,0,1]);
% xlabel('X', 'FontSize', 20);
% ylabel('Y', 'FontSize', 20);
% set(gca, 'FontSize', 20)


%%
if 1

x1 = fid1(:, 2);
y1 = fid1(:, 3);

x1_plot = x1(find(y1==0.5));
u1_ploty = u1(find(y1==0.5));
v1_ploty = v1(find(y1==0.5));
T1_ploty = temp1(find(y1==0.5));

y1_plot = y1(find(x1==0.5));
u1_plotx = u1(find(x1==0.5));
v1_plotx = v1(find(x1==0.5));
T1_plotx = temp1(find(x1==0.5));

%%%%%%%%%%%%%%%%%%
x2 = fid2(:, 2);
y2 = fid2(:, 3);

x2_plot = x2(find(y2==0.5));
u2_ploty = u2(find(y2==0.5));
v2_ploty = v2(find(y2==0.5));
T2_ploty = temp2(find(y2==0.5));

y2_plot = y2(find(x2==0.5));
u2_plotx = u2(find(x2==0.5));
v2_plotx = v2(find(x2==0.5));
T2_plotx = temp2(find(x2==0.5));

%%%%%%%%%%%%%%%%%%%
x3 = fid3(:, 2);
y3 = fid3(:, 3);

x3_plot = x3(find(y3==0.5));
u3_ploty = u3(find(y3==0.5));
v3_ploty = v3(find(y3==0.5));
T3_ploty = temp3(find(y3==0.5));

y3_plot = y3(find(x3==0.5));
u3_plotx = u3(find(x3==0.5));
v3_plotx = v3(find(x3==0.5));
T3_plotx = temp3(find(x3==0.5));

%%%%%%%%%%%%%%%%%%%
figure
plot(x1_plot, v1_ploty, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(x2_plot, v2_ploty, '--sk', 'LineWidth', 2, 'MarkerSize', 10);
plot(x3_plot, v3_ploty, '-.*k', 'LineWidth', 2, 'MarkerSize', 10); 
hold off;
xlabel('x^{*}', 'FontSize', 20);
ylabel('v^{*}', 'FontSize', 20);
legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Ra = 10^{5}');
set(gca, 'FontSize', 20)

figure
plot(x1_plot, u1_ploty, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(x2_plot, u2_ploty, '--sk', 'LineWidth', 2, 'MarkerSize', 10);
plot(x3_plot, u3_ploty, '-.*k', 'LineWidth', 2, 'MarkerSize', 10); 
hold off;
xlabel('x^{*}', 'FontSize', 20);
ylabel('u^{*}', 'FontSize', 20);
legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Ra = 10^{5}');
set(gca, 'FontSize', 20)



figure
plot(u1_plotx, y1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(u2_plotx, y2_plot, '--sk', 'LineWidth', 2, 'MarkerSize', 10)
plot(u3_plotx, y3_plot, '-.*k', 'LineWidth', 2, 'MarkerSize', 10); 
hold off
xlabel('u^{*}', 'FontSize', 20);
ylabel('y^{*}', 'FontSize', 20);
legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Ra = 10^{5}');
set(gca, 'FontSize', 20)

figure
plot(v1_plotx, y1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
plot(v2_plotx, y2_plot, '--sk', 'LineWidth', 2, 'MarkerSize', 10)
plot(v3_plotx, y3_plot, '-.*k', 'LineWidth', 2, 'MarkerSize', 10); 
hold off
xlabel('v^{*}', 'FontSize', 20);
ylabel('y^{*}', 'FontSize', 20);
legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Ra = 10^{5}');
set(gca, 'FontSize', 20)

%%%%% T  %%%%%%%%%%%%%%
% figure
% plot(x1_plot, T1_ploty, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
% plot(x2_plot, T2_ploty, '--sk', 'LineWidth', 2, 'MarkerSize', 10)
% plot(x3_plot, T3_ploty, '-.*k', 'LineWidth', 2, 'MarkerSize', 10)
% hold off
% xlabel('x^{*}', 'FontSize', 20);
% ylabel('T^{*}', 'FontSize', 20);
% axis([0,1,0,0.6]);
% legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Ra = 10^{5}');
% set(gca, 'FontSize', 20)
% 
% figure
% plot(T1_plotx, y1_plot, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on
% plot(T2_plotx, y2_plot, '--sk', 'LineWidth', 2, 'MarkerSize', 10)
% plot(T3_plotx, y3_plot, '-.*k', 'LineWidth', 2, 'MarkerSize', 10)
% hold off
% xlabel('T^{*}', 'FontSize', 20);
% ylabel('y^{*}', 'FontSize', 20);
% axis([0,1,0,1]);
% legend('Ra = 10^{3}', 'Ra = 10^{4}', 'Ra = 10^{5}');
% set(gca, 'FontSize', 20)


end





