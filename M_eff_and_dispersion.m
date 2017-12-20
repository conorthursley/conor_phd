m1=0.1; 
m2=1*m1;
k1=1000;
k2=1*k1; %118.4353;
SpringRatio=k2/k1
MassRatio=m2/m1
% m1=0.2927;
% m2=0.03;
% k1=8557;
% k2=1895;

w0=sqrt(k2/m2);

% effective mass
theta=m2/m1;
w2=sqrt(k2/m2);

w=linspace(0,1000,80000);
A=w;
B=w2;
C=A./B;

meff=m1+(m2*w2^2)./(w2^2-A.^2);

figure
plot(meff/m1,C)
grid
% axis([-20 20 0 3])

qL=2*(asin(sqrt((meff)/(4*k1).*A.^2)));

figure
plot(real(qL)/pi,A/(2*pi))
grid
% axis([0 1 0 12])
% omega=(1/B^2)*sqrt((2*k1)./((meff)).*(1-cos(qL)));
% 
% figure
% plot(qL,abs(omega),qL,imag(omega),qL,real(omega))
% grid
% axis([0 pi 0 3])

% syms wL qL %m1 m2 k1 k2
% f1= @(qL,wL) 0.03*wL^4 - (m1+m2)*k2 + 2*m2*k1*(1-cos(qL))*wL^2 + 2*k1*k2*(1-cos(qL));
% % sol = solve(eqn)
% % solw=solve(eqn,wL)
% % a=linspace(0,pi,10000);
% figure
% fimplicit(f,[0 pi 0 3.5])
% xlabel('wavenumber, qL') % x-axis label
% ylabel('frequency, w') % y-axis label
% grid
% omega=sqrt(((4*k1)/m1))*sin(qL/2);
% Q1=-(10000*cos(a) - (50*(360000*cos(a).^2 - 717600.*cos(a) + 8940012/25).^(1/2))/3 - 10000);
% Q2=-(10000*cos(a) + (50*(360000*cos(a).^2 - 717600.*cos(a) + 8940012/25).^(1/2))/3 - 10000);
% 
% figure
% plot(a,(Q1/50),a,(Q2/50))
% grid
% % axis([0 pi 0 3])
% 
%% Transmittance

% T10=(k1)./(k1-meff.*A.^2);
% figure
% loglog(A,T10)
% grid
% hold on
% T9=(k1)./(k1*(2-T10)-meff.*A.^2);
% loglog(A,T9)
% T8=(k1)./(k1*(2-T9)-meff.*A.^2);
% loglog(A,T8)
figure
hold on
grid
B=A/(2*pi);
T7=(k1)./(k1*(1)-meff.*A.^2);
loglog(B,T7)
T6=(k1)./(k1*(2-T7)-meff.*A.^2);
loglog(B,T6)
T5=(k1)./(k1*(2-T6)-meff.*A.^2);
loglog(B,T5)
T4=(k1)./(k1*(2-T5)-meff.*A.^2);
loglog(B,T4)
T3=(k1)./(k1*(2-T4)-meff.*A.^2);
loglog(B,T3)
T2=(k1)./(k1*(2-T3)-meff.*A.^2);
loglog(B,T2)
T1=(k1)./(k1*(2-T2)-meff.*A.^2);
loglog(B,T1)
figure
Tall=T7+T6+T5+T4+T3+T2+T1;
loglog(abs(Tall),B)
grid
%% band structure
freq=50;
w_2=freq*2*pi;
m2=0.03; m2
k2=w_2^2*m2; k2
upper=(freq+2)/freq;
lower=(freq-2)/freq;
Mr=upper^2-1;
m1=m2/Mr; m1

L=12500000*lower^2*(-49999041*Mr+50000000*lower^2-49999041);
Q=49999041*(50000000*lower^2-49999041);

T=L/Q;

Kr=Mr/T;
k1=k2/Kr; k1