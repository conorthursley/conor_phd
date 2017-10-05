%% Dispersion.m
% Purpose - solve for the dispersion properties of a metamaterial system
% with k1,k2,m1,m2 and a driving frequency
%--------------------------------
% Analytical equation
%----------------------------
% cosh(qa)=1+0.5*(w^4-(1-B)*K^2*w^2)/(K^2-w^2)
%-------------------------------
% qa = wavenumber x unit spacing
% B = mass ratio
% K = sqrt(A/B)
% A = spring ratio
k1=2e6;
m1=1;
k2=0.1*k1;
m2=9*m1;

B=m2/m1;
A=k2/k1;
K=sqrt(A/B);
w=0:1:50;


K=1;
B=0.1;
A=0.1;
% for i=1:length(w)
%     Q=acosh(1+0.5*((w(i)^4-(1-B)*K^2*w(i)^2)/(K^2-w(i)^2)));
%     X(i)=real(Q);
%     Y(i)=imag(Q);
% end

% plot(X,w/(2*pi),Y,w/(2*pi))
Q=linspace(0,6,100);
w1 =-(1/sqrt(2))*sqrt(sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
w2=(1/sqrt(2))*sqrt(sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
w3=-(1/sqrt(2))*sqrt(-sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
w4=(1/sqrt(2))*sqrt(-sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
w3=imag(w3);
w4=imag(w4);
plot(Q,w1,Q,w2,Q,w3,Q,w4)


grid on
xlabel('wavenumber, qL') % x-axis label
ylabel('frequency, Hz') % y-axis label
title('kill me')