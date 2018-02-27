%syms m1 m2 m3 m4 z1 z2 z3 z4 k1 k2 f1 f2 f3 f4
% A=[s^2*m1+k1+k2 -k2 -k1 0;-k1 m2*s^2+k2 0 0;-k1 0 m3*s^2+k2+k2 -k2;0 0 -k2 m4*s^2+k2];
% F=[f1;f2;f3;f4];
% Z=[z1;z2;z3;z4];
% X=[z1^2;z2^2;z3^2];

%% parameters of the AMM system
m1=1;
m2=0.3*m1;
m3=m1;
m4=m2;
k1=1000;
k2=0.1*k1;
%--------------------------------------------------------------
syms s z1 z2 z3 z4 f1 f2 f3 f4
A=[s^2*m1+k1+k2 -k2 -k1 0;-k1 m2*s^2+k2 0 0;-k1 0 m3*s^2+k2+k2 -k2;0 0 -k2 m4*s^2+k2];
F=[f1;f2;f3;f4];
Z=[z1;z2;z3;z4];
[i h u t]=solve(A*Z==F);
[n,d]=numden(i);

% then set all other expressions (f2,f3) to 0 and divide everything by f1
n=subs(n,f2,0); % then f3 and f4