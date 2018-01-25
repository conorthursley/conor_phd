%% test
m1=0.1; 
m2=0.5*m1;
k1=1000;
k2=1.5*k1;
c2=0.02;
c1=c2;

% numerator
num=[m2 c2 k2];

% denominator
den=[m1*m2 m2*c1+m1*c2 m1*k2+2*c1*c2+m2*(k1+k2) k2*c1+c2*(k1+k2) k2*(k1+k2)-k2^2];

% TF
H=tf(num,den)

% State Space form 
[A,B,C,D]=tf2ss(num,den)

Fs = 5;
dt = 1/Fs;
N = 50;
t = dt*(0:N-1);

ux = [1 zeros(1,N-1)];
u0 = zeros(1,N);
u = [ux;u0];

x = [0;0;0;0];
for k = 1:N
    y(:,k) = C*x + D*u(:,k);
    x = A*x + B*u(:,k)';
end

stem(t,y','.')
xlabel('t')
legend('a_1','a_2')
title('Mass 1 Excited')
grid
