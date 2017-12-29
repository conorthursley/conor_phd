%% force displacement curve for nonlinear stiffness
function F = piecewiseCurve(k1,k2)
kn=k2;
k1=k1;
dy=1;
y=(-10:dy:10)';  %column matrix
f=k1*y;  %linear case


syms q
spring=piecewise(q<=-5,k1,-5<q<5, k1*5*q+kn*q.^3, q>=5, k1);
k=subs(spring, q, y);
k3=double(k);


% f1=+k1*y+kn*y.^3;  %nonlinear Hardening
% f2=k1*y-kn*y.^3;  %nonlinear softening
F=[y k3];
figure(1)
plot(y,f,y,k3)
grid on
xlabel('y(m)','FontSize', 14)
ylabel('Force(N)','FontSize', 14)
title ('Force displacement curve for nonlinear function','FontSize',18)
legend({'linear case','nonlinear Hardening','nonlinear softening'},'FontSize', 14)
end

