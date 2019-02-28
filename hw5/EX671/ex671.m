clc;clear all;close all;
%%
if 0
TL88CN = load('TL88CN.dat');
TL88BD = load('TL88BD.dat');
TL88FD = load('TL88FD.dat');

TQ44CN = load('TQ44CN.dat');
TQ44BD = load('TQ44BD.dat');
TQ44FD = load('TQ44FD.dat');

time = TL88CN(:,1);

T1CN = TL88CN(:,2);
T1BD = TL88BD(:,2);
T1FD = TL88FD(:,2);

T2CN = TQ44CN(:,2);
T2BD = TQ44BD(:,2);
T2FD = TQ44FD(:,2);

figure
plot(time, T1CN, '-ok', 'LineWidth', 2, 'MarkerSize', 10); hold on; grid on
plot(time, T1BD, '--sk', 'LineWidth', 2, 'MarkerSize', 10);
plot(time, T1FD, '-.+k', 'LineWidth', 2, 'MarkerSize', 15); hold off
axis([0,1.4, 0,0.35])
set(gca, 'xtick',0:0.2:1.4, 'ytick', 0:0.05:0.35);

xlabel('Time, t', 'FontSize', 20)
ylabel('Temperature, {\theta}(0,0,t)', 'FontSize', 20)
legend('Crank-Nicolson scheme', 'Backward difference', 'Forward difference','Location','SouthEast')
set(gca, 'FontSize', 20)
end

%%
if 1
fcn = load('ex671CN.txt');
fbd = load('ex671BD.txt');

x = fcn(1:9, 2);
T_t025_cn = fcn(1:9, end);
T_t050_cn = fcn(10:18, end);
T_t100_cn = fcn(19:27, end);
T_t125_cn = fcn(28:36, end);

T_t025_bd = fbd(1:9, end);
T_t050_bd = fbd(10:18, end);
T_t100_bd = fbd(19:27, end);
T_t125_bd = fbd(28:36, end);

figure
plot(x, T_t025_cn, '-ok', 'LineWIdth', 2, 'MarkerSize', 10);grid on; hold on
plot(x, T_t050_cn, '-sk', 'LineWIdth', 2, 'MarkerSize', 10);
plot(x, T_t100_cn, '-*k', 'LineWIdth', 2, 'MarkerSize', 10);
plot(x, T_t125_cn, '-+k', 'LineWIdth', 2, 'MarkerSize', 10);

plot(x, T_t025_bd, '--ok', 'LineWIdth', 2, 'MarkerSize', 10);
plot(x, T_t050_bd, '--sk', 'LineWIdth', 2, 'MarkerSize', 10);
plot(x, T_t100_bd, '--*k', 'LineWIdth', 2, 'MarkerSize', 10);
plot(x, T_t125_bd, '--+k', 'LineWIdth', 2, 'MarkerSize', 10);hold off;

axis([0,1.0, 0,0.35])
set(gca, 'xtick',0:0.1:1.0, 'ytick', 0:0.05:0.35);
xlabel('Distance, x', 'FontSize', 20)
ylabel('Temperature, {\theta}(x,0,t)', 'FontSize', 20)
legend('{\alpha}=0.5, t=0.25', '{\alpha}=0.5, t=0.50', '{\alpha}=0.5, t=1.00', '{\alpha}=0.5, t=1.25',...
    '{\alpha}=1.0, t=0.25', '{\alpha}=1.0, t=0.50', '{\alpha}=1.0, t=1.00', '{\alpha}=1.0, t=1.25');
set(gca, 'FontSize', 20)

end




