clc;clear all;close all;

%%
if 0
f1 = load('EX1081Q53_U.txt');
f2 = load('EX1081Q106_U.txt');

x1 = f1(:, 2);
y1 = f1(:, 3);
u1 = f1(:, end-1);
y1_plot1 = y1(find(x1==4));
y1_plot2 = y1(find(x1==6));

u1_4 = u1(find(x1==4));
u1_6 = u1(find(x1==6));

x2 = f2(:,2);
y2 = f2(:,3);
u2 = f2(:, end-1);
y2_plot1 = y2(find(x2==4));
y2_plot2 = y2(find(x2==6));

u2_4 = u2(find(x2==4));
u2_6 = u2(find(x2==6));


% exact soln
Y = 0:0.1:2;
X1 = 4;
X2 = 6;
u_exact1 = 3/4*X1*(1-Y.^2/4);
u_exact2 = 3/4*X2*(1-Y.^2/4);

figure
plot(u1_4, y1_plot1, 'ok', 'LineWidth', 3, 'MarkerSize', 10, 'markerfacecolor', [0,0,0]);hold on; grid on
plot(u1_6, y1_plot2, 'sk', 'LineWidth', 3, 'MarkerSize', 10, 'markerfacecolor', [0,0,0]); 
plot(u2_4, y2_plot1, 'ok', 'LineWidth', 3, 'MarkerSize', 10);
plot(u2_6, y2_plot2, 'sk', 'LineWidth', 3, 'MarkerSize', 10);

plot(u_exact1, Y, '-k', 'LineWidth', 2);
plot(u_exact2, Y, '--k', 'LineWidth', 2); hold off

xlabel('Horizontal velocity, u(x,y)', 'FontSize', 20);
ylabel('Distance, y', 'FontSize', 20);
legend('x = 4, 5x3Q9', 'x = 6, 5x3Q9', 'x = 4, 10x6Q9', 'x = 6, 10x6Q9', 'x = 4, exact', 'x = 6, exact');

axis([0,4.5,0,2]);
set(gca,'xtick',0:0.5:4.5,'ytick',0:0.5:2)
set(gca, 'FontSize', 20);

end

%%
if 1
f1 = load('EX1081Q53_P.txt');
f2 = load('EX1081Q106_P.txt');

x1 = f1(:,1);
y1 = f1(:,2);
pres1 = f1(:,end);

y_top1 = max(y1);
y_bot1 = min(y1);

x1_plot1 = x1(find(y1==y_top1));
pres1_plot1 = pres1(find(y1==y_top1));

x1_plot2 = x1(find(y1==y_bot1));
pres1_plot2 = pres1(find(y1==y_bot1));

x2 = f2(:,1);
y2 = f2(:,2);
pres2 = f2(:,end);

y_top2 = max(y2);
y_bot2 = min(y2);


x2_plot1 = x2(find(y2==y_top2));
pres2_plot1 = pres2(find(y2==y_top2));

x2_plot2 = x2(find(y2==y_bot2));
pres2_plot2 = pres2(find(y2==y_bot2));

% exact soln
X1 = 0:0.1:6;
Y_top = 2;
Y_bot = 0;
pres1_exact = 3/16*(6^2 + Y_top^2 - X1.^2);
pres2_exact = 3/16*(6^2 + Y_bot^2 - X1.^2);


figure
plot(x1_plot1, pres1_plot1, '-ok', 'LineWidth', 2, 'MarkerSize', 10);hold on; grid on
plot(x2_plot1, pres2_plot1, '--sk','LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0,0,0])

plot(X1, pres1_exact, '-k', 'LineWidth', 2); hold off
xlabel('Distance, x', 'FontSize', 20);
ylabel('Pressure, P', 'FontSize', 20);
legend('5x3Q9', '10x6Q9', 'Exact solution');
axis([0,6,0,13]);
set(gca,'xtick',0:1:6,'ytick',0:1:13)
set(gca, 'FontSize', 20);

figure
plot(x1_plot2, pres1_plot2, '-ok', 'LineWidth', 2, 'MarkerSize', 10);hold on; grid on
plot(x2_plot2, pres2_plot2, '--sk', 'LineWidth', 2, 'MarkerSize', 10, 'markerfacecolor', [0,0,0]);

plot(X1, pres2_exact, '-k', 'LineWidth', 2); hold off
xlabel('Distance, x', 'FontSize', 20);
ylabel('Pressure, P', 'FontSize', 20);
legend('5x3Q9', '10x6Q9',  'Exact solution');
axis([0,6,-1,8]);
set(gca,'xtick',0:1:6,'ytick',-1:1:8)
set(gca, 'FontSize', 20);

end





