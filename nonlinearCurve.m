%% force displacement curve for nonlinear stiffness
function F = nonlinearCurve(k1,k2)
kn=100*1500;
k1=1500;
dy=0.1;
y=(-1:dy:1)';  %column matrix
f=k1*y;  %linear case
f1=+k1*y+kn*y.^3;  %nonlinear Hardening
f2=k1*y-kn*y.^3;  %nonlinear softening
F=[y f1];
figure(1)
plot(y,f,y,f1) %,y,f2)
grid on
xlabel('y(m)','FontSize', 14)
ylabel('Force(N)','FontSize', 14)
title ('Force displacement curve for nonlinear function','FontSize',18)
legend({'linear case','nonlinear Hardening','nonlinear softening'},'FontSize', 14)
end

