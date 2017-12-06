%% force displacement curve for nonlinear stiffness

clear;
k2=1e3;
kn=0.01*k2;
dy=1;
y=(-10:dy:10)';  %column matrix
f=k2*y;  %linear case
f1=k2*y+kn*y.^3;  %nonlinear Hardening
f2=k2*y-kn*y.^3;  %nonlinear softening
F=[y f2];
figure(1)
plot(y,f,y,f1,y,f2)
grid on
xlabel('y(m)','FontSize', 14)
ylabel('Force(N)','FontSize', 14)
title ('Force displacement curve for nonlinear function','FontSize',18)
legend({'linear case','nonlinear Hardening','nonlinear softening'},'FontSize', 14)

