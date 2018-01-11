%% Piecewise Curve
function F=piecewiseCurve(k1,k2,s)
k2; %middle of the restoring force curve
k1; %edges of the restoring force curve
% when K2 is 0, the curve should be linear
% the larger ratio (k1/k2) is, the more sharper the FD curve is
% our distance function
dy=0.01;
y=(-0.1:dy:0.1)';  %column matrix
f=k2*y;  %linear case
f1=(y<-s).*((k1)*y+(k1-k2)*s)+(y>s).*(k1*y+(k2-k1)*s)+(y<=s).*(y>=-s).*(k2*y);
F=[y f1];
figure(1)
plot(y,f,y,f1)
grid on
xlabel('y(m)','FontSize', 14)
ylabel('Force(N)','FontSize', 14)
title ('Force displacement curve for nonlinear function','FontSize',18)
legend({'linear case','Piecewise'},'FontSize', 14)