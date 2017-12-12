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

figure
plot(C,meff)
grid

qL=acos(1-(abs(meff)./(2*k1)).*A.^2);

figure
plot(real(qL),C,imag(qL),C)
grid
