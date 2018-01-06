%% Piecewise Curve
k2=240; %middle of the restoring force curve
k1=1000; %edges of the restoring force curve
s=4; % our distance function
dy=1;
y=(-10:dy:10)';  %column matrix
f=k1*y;  %linear case
f1=(y<-s).*((k1)*y+(k1-k2)*s)+(y>s).*(k1*y+(k2-k1)*s)+(y<=s).*(y>=-s).*(k2*y);
F=[y f1];
figure(1)
plot(y,f,y,f1)
grid on
xlabel('y(m)','FontSize', 14)
ylabel('Force(N)','FontSize', 14)
title ('Force displacement curve for nonlinear function','FontSize',18)
legend({'linear case','Piecewise'},'FontSize', 14)