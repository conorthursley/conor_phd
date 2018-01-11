%% Duffing Oscillator solved using ode45
% use ode45 to numerically solve the duffing oscillator ODE equation
% x2 + delta*x1 + beta*x0 + alpha*x0^3 = gamma*cos(omega*t)
% x2= double derivative of x wrt to t.
%----------------------------------------------

%% Coefficients 
alpha=1;
beta=5;
delta=0.02;
gamma=8;
omega=0.5;
tspan=[0 100];
y0=[0 0]; %initial conditions
%% ODE solver
tic
[t, result]=ode45(@(t,y) duffing(t, y, alpha, beta, gamma, delta, omega), tspan, y0);
toc


%% Plots
figure;
plot(t,result);
xlabel('Time');
ylabel('State');
grid on
a=num2str(alpha);
str ='Forced duffing oscillator of x^{''}+\deltax1 + \alphax0 + \betax0^3 = \gammacos(wt) with \alpha=a, \beta=%d, \gamma=%d, \delta=%d, \omega=%d',alpha,beta,gamma,delta,omega;
title(str,'FontSize',12)
legend({'disp','velo'},'FontSize',14)


figure;
plot(result(:,1),result(:,2));
title(str,'FontSize',14)
xlabel('Displacement');
ylabel('Velocity');
grid on


%% Function

function dy = duffing(t, y, alpha, beta, gamma, delta, omega)
A=y(2);
B=gamma*cos(omega*t)-delta*y(2)-beta*y(1)-alpha*y(1).^3;
dy=[A; B];
end
