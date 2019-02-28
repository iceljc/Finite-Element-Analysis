clc;clear all;

% f1 = load('RE0_L88.txt');
% f2 = load('RE250_L88.txt');
% f3 = load('RE500_L88.txt');
% f4 = load('RE750_L88.txt');

f1 = load('RE0_Q44.txt');
f2 = load('RE250_Q44.txt');
f3 = load('RE500_Q44.txt');
f4 = load('RE750_Q44.txt');

x = f1(:,2);
y = f1(:,3);

u1 = f1(:,4);
u2 = f2(:,4);
u3 = f3(:,4);
u4 = f4(:,4);

y_table = y(find(x==0.5));
u1_table = u1(find(x==0.5));
u2_table = u2(find(x==0.5));
u3_table = u3(find(x==0.5));
u4_table = u4(find(x==0.5));

u_L88 = [y_table, u1_table, u2_table, u3_table, u4_table];
u_L88 = u_L88(2:end-1,:);