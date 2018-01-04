tspan = [0 10];
k = 1e3;
c = 40;
k3 = 1.1e5;
A = 30;
omega = 10;

ode_fun = @ (t,u,k,c,k3,A,omega
[u(2);-k*u(1)-c*u(2)-k3*u(1).^3+A*sin(omega*t)];

[t,u] = ode15s(@(t,u)ode_fun(t,u,k,c,k3,A,omega),tspan,[0 0]);

figure;
plot(t,u);
xlabel('Time');
ylabel('State');

acc = -k*u(:,1)-c*u(:,2)-k3*u(:,1).^3+A*sin(omega*t);
figure;
plot(u(:,1),acc);
xlabel('Position');
ylabel('Acceleration');