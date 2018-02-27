%% ODE computation of single mass spring system - BASIC


tspan = [0 1000];
y0=[0,0];

[t,y]=ode45(@fun, tspan, y0);
figure
plot(t,y(:,1))

function dy=fun(t,y)
x1=y(1);
x2=y(2);
m=0.1;
k=1000;
F=0.1*cos(0.01*t);
dx1=x2;
dx2=(F-2*k*x1)/m;
dy=[dx1;dx2];
end


    