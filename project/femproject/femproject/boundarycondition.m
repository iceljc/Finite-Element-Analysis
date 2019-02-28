clc;clear all;

a = 1:1:101;     a = a';
b = 101:101:10101;  b = b';
c = 202:101:10201;  c = c';
d = 10102:1:10200; d = d';

%%
A = [a , ones(size(a,1), 1)];
A1 = reshape(A', 1, 2*size(a,1));

A = [a , 2*ones(size(a,1), 1)];
A2 = reshape(A', 1, 2*size(a,1));

%%
B = [b , ones(size(b,1), 1)];
B1 = reshape(B', 1, 2*size(b,1));

B = [b , 2*ones(size(b,1), 1)];
B2 = reshape(B', 1, 2*size(b,1));

%% 
C = [c , ones(size(c,1), 1)];
C1 = reshape(C', 1, 2*size(c,1));

C = [c , 2*ones(size(c,1), 1)];
C2 = reshape(C', 1, 2*size(c,1));

%% 
D = [d , ones(size(d,1), 1)];
D1 = reshape(D', 1, 2*size(d,1));

D = [d , 2*ones(size(d,1), 1)];
D2 = reshape(D', 1, 2*size(d,1));

%%

bc1 = [A1, A2, B1, B2, C1, C2, D1, D2];

BC1 = reshape(bc1, 16, 100)';

%%

temp1 = 1:101:10101;    temp1 = temp1';
temp2 = 101:101:10201;  temp2 = temp2';

T = [temp1 , 3*ones(size(temp1,1), 1)];
T1 = reshape(T', 1, 2*size(temp1,1));

T = [temp2 , 3*ones(size(temp2,1), 1)];
T2 = reshape(T', 1, 2*size(temp2,1));

bc2 = [T1, T2, zeros(1,12)];
BC2 = reshape(bc2, 16, 26)';

vspv = [ones(1, 101), zeros(1, 101), zeros(1,6)];

VSPV = reshape(vspv, 8, 26)';









