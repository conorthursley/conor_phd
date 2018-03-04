%% ODE computation of single mass spring system - BASIC

clear all
tspan = [0 12];
y0=[5e-3,(2*sqrt(2)*10^(-3))];

[t,y]=ode45(@fun, tspan, y0);
figure
plot(t,y)

function dy=fun(t,y)
x1=y(1);
x2=y(2);
% m=0.1;
% k=1000;
F=1; %*cos(0.01*t);
dx1=x2;
dx2=(1+(-4)*x1);
dy=[dx1;dx2];
end


    