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

qL=linspace(0,pi,10000);
% omega=(1/B^2)*sqrt((2*k1)./((meff)).*(1-cos(qL)));
% 
% figure
% plot(qL,abs(omega),qL,imag(omega),qL,real(omega))
% grid
% axis([0 pi 0 3])

syms well qLel m k
eqn = cos(qLel) == 1-(m/(2*k))*well^2;
sol = solve(eqn)
solw=solve(eqn,well)

omega=(2^(1/2)*(-k1)^(1/2)*(cos(qL) - 1).^(1/2))./abs(meff).^(1/2);

figure
plot(qL,abs(omega),qL,imag(omega),qL,real(omega))
grid
axis([0 pi 0 3])

