function [dy]=piecewise(t,y,input,k1,m1,k2,m2)

% Define parameters of spring mass model
% k1=2e6;
% m1=1;
% k2=0.1*k1;
% m2=0.05*m1;


% define X matrix, Xdot
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);

syms q
spring=piecewise(q<=5e-4,k1,5e-4<q<1.5e-3, k2, q>=1.5e-3, k1);
k=subs(spring, q, y1);
k3=double(k);
% k3=k2;
% w0=sqrt(k2/m2);

f=input(1)*(2*pi); %Hz to rad/s
% Define A matrix
A=[x2; ...
    -((k1/m1)+(k3/m1))*(x1)+((k3/m1)*y1); ...
    y2; ...
    ((k3/m2)*x1)-((k3/m2)*y1)];

% Input excitation force
%sinusodial harmonic
H=input(2)*(sin(f*t));

%sawtooth with random noise
% s=sawtooth(t);
% noise=awgn(s,0.5);   
% 
% % random noise
% R=0.1*randn(size(t));
B = [0 H/m1 0 0];


%output result as the equation 
dy=A+B';
% dy=A*y;

end






