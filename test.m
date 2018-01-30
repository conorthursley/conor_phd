%% test
m1=0.1; 
m2=0.5*m1;
k1=1000;
k2=1.5*k1;
c2=0.00;
c1=c2;

% numerator
num=[m2 c2 k2];
% num=1;

% denominator
den=[m1*m2 m2*c1+m1*c2 m1*k2+2*c1*c2+m2*(k1+k2) k2*c1+c2*(k1+k2) k2*(k1+k2)-k2^2];

% TF
H=tf(num,den);
options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
figure(1)
bode(H,options);
grid on
set(findall(gcf,'Type','text'),'FontSize',20)
set(findall(gcf,'Type','axes'),'FontSize',20,'LineWidth',1,'XColor','black','YColor','black')

%% 2 cell system
% numerator
num2=[m2^2 0 2*m2*k2 0 k2^2];
% denominator 
den2=[m1*m1*m2*m2 0 2*m1*m2*m2*(k1+k2) 0 2*m1*m2*k2*k1+2*m2*m2*k1*k2+m2*m2*k1*k1+m2*m2*k2*k2 0 2*m2*(k1+k2)*k2*k1-k1*k1*m2 0 k2*k2*k1*k1-k1*k1*k2];

H1=tf(num2,den2)

options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
figure(1)
bode(H1,options);
grid on
set(findall(gcf,'Type','text'),'FontSize',20)
set(findall(gcf,'Type','axes'),'FontSize',20,'LineWidth',1,'XColor','black','YColor','black')
%% SS
[A,B,C,D]=tf2ss(num,den)

Fs = 5;
dt = 1/Fs;
N = 50;
t = dt*(0:N-1);

ux = [1 zeros(1,N-1)];
u0 = zeros(1,N);
u = [ux;u0];

x = [0;0;0;0];
for k = 1:N
    y(:,k) = C*x + D*u(:,k);
    x = A*x + B*u(:,k)';
end
[b,a] = ss2tf(A,B,C,D);
G=H(b,a)

bode(G)


