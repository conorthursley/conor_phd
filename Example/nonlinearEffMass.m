clear all
close all
m=1;
kl=1e+3; %linear stiffness, k2
kn=1; %nonlinear stiffness
c=1; %damping
f=200; %forcing amplitude 
m1=0.1;
m2=0.9;
% w=150;
mst=m1+m2;
%-------------------------------------------------
w0=sqrt(kl/m2);
A=10;

f = @(omega,Y) sqrt(((64/(243*(kn^3)))*(((kl^3)/3)-((m2^3*w0^6*omega^6)/3)-(kl^2*m2*w0^2*omega^2)+(kl*m2^2*w0^4*omega^4)))+(4*m2^2*w0^4*omega^4*A^2)/(9*(kn^2)))+(2*m2*w0^2*omega^2*A)/(3*kn)-(Y^3);

figure
fimplicit(f,[0 3.5 -100 100])
xlabel('wavenumber, qL') % x-axis label
ylabel('frequency, w') % y-axis label
title('Huang et al. 2009')
