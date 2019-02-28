clc;clear all;

syms x y
f1 = 2/3*x*(x + 1/2)*(x - 1/2)*(x - 1);
f2 = -8/3*(x + 1)*x*(x - 1/2)*(x - 1);
f3 = 4*(x + 1)*(x + 1/2)*(x - 1/2)*(x - 1);
f4 = -8/3*(x + 1)*(x + 1/2)*x*(x - 1);
f5 = 2/3*(x + 1)*(x + 1/2)*x*(x - 1/2);

g1 = 2/3*y*(y + 1/2)*(y - 1/2)*(y - 1);
g2 = -8/3*(y + 1)*y*(y - 1/2)*(y - 1);
g3 = 4*(y + 1)*(y + 1/2)*(y - 1/2)*(y - 1);
g4 = -8/3*(y + 1)*(y + 1/2)*y*(y - 1);
g5 = 2/3*(y + 1)*(y + 1/2)*y*(y - 1/2);

F1 = [f1; f2; f3; f4; f5];
F2 = [g1, g2, g3, g4, g5];

phi = F1*F2;

Phi = simplify(expand(phi));

DSFL1 = diff(phi, x);
DSFL2 = diff(phi, y);

DSFL1 = simplify(DSFL1);
DSFL2 = simplify(DSFL2);





