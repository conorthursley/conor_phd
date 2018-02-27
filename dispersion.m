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
%% First attempt

% B=m2/m1;
% A=k2/k1;
% K=sqrt(A/B);
% w=0:1:50;

% for i=1:length(w)
%     Q=acosh(1+0.5*((w(i)^4-(1-B)*K^2*w(i)^2)/(K^2-w(i)^2)));
%     X(i)=real(Q);
%     Y(i)=imag(Q);
% end

% % plot(X,w/(2*pi),Y,w/(2*pi))
% A=linspace(0,2*pi,100);
% w1 =-(1/sqrt(2))*sqrt(sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
% w2=(1/sqrt(2))*sqrt(sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
% w3=-(1/sqrt(2))*sqrt(-sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
% w4=(1/sqrt(2))*sqrt(-sqrt((B*K^2 - K^2 + 2*cosh(Q) - 2).^2 - 4*(2*K^2 - 2*K^2*cosh(Q))) - B*K^2 + K^2 - 2*cosh(Q) + 2);
% w3=imag(w3);
% w4=imag(w4);
% plot(Q,w1,Q,w2,Q,w3,Q,w4)


% grid on
% xlabel('wavenumber, qL') % x-axis label
% ylabel('frequency, Hz') % y-axis label
% title('kill me')
%% Second Attempt
% for i=1:length(A)
%     w1(i)=+(sqrt(-1)*sqrt(k2)*sqrt(2*k1*cos(A(i))-2*k1+m1))/(sqrt(m2*(-2*k1*cos(A(i))+2*k1+k2+m1)));
%     w2(i)=-w1(i);
% end
% 
% plot(A,w1,'r+',A,w2,'bo')

%% Third Attempt
% w=linspace(0,pi,100);
% for i = 1:length(w)
%     qL(i,:) = acos(1-((w(i)^2/2*k1))*m1-((k2*m2)/(m2*w(i)^2-k2)));
% 
% end
% 
% X=real(qL);
% Y=imag(qL);
% plot(X,w,'r',Y,w,'b')
% grid on

%% Fourth Attempt
% w is angular freqeuncy and changes
% w0 is our resonance frequency and is w0=sqrt(k^2/m^2)
w0=sqrt(k2/m2); % rad/s
% qL is the dimensionless wave number. in this case, the values will take
qL=linspace(0,pi,100); %evenly spaced vector from 0 to pi
A=qL;
% creat w12 which means roots 1 and 2 of angular frequency, w.
% this value takes on a + and - value, hence 1 and 2
% this value was created using Wolfram Alpha computational solver for
% quartics with the equation used as f in the above section.
w12=sqrt(-sqrt((-2*m1*k1*(1-cos(A))-k2*m1-k2*m2).^2-4*m1*m2*(2*k1*k2-2*k1*k2*cos(A)))/(m1*m2)+(2*k1*(1-cos(A)))/m1+k2/m1+k2/m2)/sqrt(2);
% need to make the frequency ratio for the dispersion curve
Omega1=(w12)/w0;

% creat w34 which means roots 3 and 4 of angular frequency, w.
% this value takes on a + and - value, hence 3 and 4
w34=-sqrt(sqrt((-2*m1*k1*(1-cos(A))-k2*m1-k2*m2).^2-4*m1*m2*(2*k1*k2-2*k1*k2*cos(A)))/(m1*m2)+(2*k1*(1-cos(A)))/m1+k2/m1+k2/m2)/sqrt(2);
% need to make the frequency ratio for the dispersion curve
Omega2=(w34)/w0;

figure
plot(A,+Omega1,A,+Omega2)
grid on
xlabel('wavenumber, qL') % x-axis label
ylabel('frequency ratio, w/w_0') % y-axis label
title('Analytical determination of Dispersion Curve')
