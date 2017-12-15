m1=0.1; 
m2=0.3;
k1=1000;
k2=10;
w0=sqrt(k2/m2);

% effective mass
theta=m2/m1;
w2=sqrt(k2/m2);

w=linspace(0,15,10000);
A=w;
B=w2;
C=A./B;

meff=m1*(C-theta-1)./(C-1);

% figure
% plot((meff),C)
% grid
% axis([-20 20 0 3])

% qL=linspace(0,pi,10000);
% omega=(1/B^2)*sqrt((2*k1)./((meff)).*(1-cos(qL)));
% 
% figure
% plot(qL,abs(omega),qL,imag(omega),qL,real(omega))
% grid
% axis([0 pi 0 3])

% syms wL qL %m1 m2 k1 k2
f1= @(qL,wL) 0.03*wL^4 - (m1+m2)*k2 + 2*m2*k1*(1-cos(qL))*wL^2 + 2*k1*k2*(1-cos(qL));
% sol = solve(eqn)
% solw=solve(eqn,wL)
% a=linspace(0,pi,10000);
figure
fimplicit(f,[0 pi 0 3.5])
xlabel('wavenumber, qL') % x-axis label
ylabel('frequency, w') % y-axis label
grid
% omega=sqrt(((4*k1)/m1))*sin(qL/2);
% Q1=-(10000*cos(a) - (50*(360000*cos(a).^2 - 717600.*cos(a) + 8940012/25).^(1/2))/3 - 10000);
% Q2=-(10000*cos(a) + (50*(360000*cos(a).^2 - 717600.*cos(a) + 8940012/25).^(1/2))/3 - 10000);
% 
% figure
% plot(a,(Q1/50),a,(Q2/50))
% grid
% % axis([0 pi 0 3])
% 
