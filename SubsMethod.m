%% Substituion method for solving duffing oscillator
%------------------------
%Negative effective mass in acoustic metamaterial with nonlinear mass-in-mass subsystems
% L. Cveticanin 2017
%------------------------

% X=u1-u2

m1=1;
m2=5;
k1=1e5;
k2=1e4;
k3=1e3;

%time step and initial condition
tspan = 0:0.001:10;
y10 = 0; y20 = 0; y30=0; y40=0;
y0 = [y10; y20; y30; y40];
op=odeset('abstol',1e-9,'reltol',1e-9);
A=1; %amplitude of exciting force
freq=5; %Hz
w=freq*2*pi;
F=A*sin(freq*2*pi*tspan);
% (tspan,y0,)
[t,y] = ode45(@(t,y) kms(t,y0,m1,m2,k2,k3,A,w),tspan,y0,op);

x1=y(:,2); x2=y(:,4);
% plot(x1,x2);  %plot the variable x and y
figure
plot(t,x1,t,x2)
legend 'x1' 'x2'
% reconfigure X and Y to get u1 and u2 seperately
% X=x1, Y=x2
u1=x1+x2;
u2=x2;
figure
plot(t,u1,t,u2)

function dy = kms(t,y0,m1,m2,k2,k3,A,w)

x1=y0(1);
x2=y0(2);
y1=y0(3);
y2=y0(4);

dx1=x2;
dx2=A*sin(w*t)/m1 - ((m1+m2)*k2*x1)/(m1*m2)-((m1+m2)*k3*x1^3)/(m1*m2);
dx3=y2;
dx4=((A*sin(w*t))-m1*x2)/(m1+m2);
dy = [dx1; dx2; dx3; dx4];
end


% amp=0.42;  % control parameter
% b=0.5;
% alpha=-1.0d0;
% beta=1.0d0;
% w=1.0;
% %time step and initial condition
% tspan = 0:0.1:500;
% x10 = 0.5021; x20 = 0.17606;
% y0 = [x10; x20];
% op=odeset('abstol',1e-9,'reltol',1e-9);
% [t,y] = ode45(@(t,x) f(tspan,y0,b,alpha,beta,amp,w),tspan,y0,op) ;
% x1=y(:,1); x2=y(:,2);
% plot(x1,x2);  %plot the variable x and y
% fprintf(fid,'%12.6f\n',x1,x2);
% 
% function dy = f(t,y,b,alpha,beta,amp,w)
% x1 = y(1);    x2 = y(2);
% dx1=x2;
% dx2=-b*x2-alpha*x1-beta*x1^3+amp*sin(w*t);
% dy = [dx1; dx2];
% end
% 
