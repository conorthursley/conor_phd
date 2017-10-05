function Xdot=Xd(t,y)
%Define parameters from model
k1=2e6;
m1=1;
k2=2e5;
m2=9*m1;
k3=2e4;

%define 'a' coefficient to resolve large solution issue
a=1e-6;
%define matrices from initial conditions
x1=a*y(1);
x2=y(2);
y1=a*y(3);
y2=y(4);

%set nonlinear term, k3, to 0
k3=0;

%define A matrix for EOM
A=[x2; -((k2/m1)*(x1-y1))-((k3/m1)*((x1-y1)^3))+((k1/m1)*x1); y2; -((k2/m2)*(y1-x1))-((k3/m2)*((y1-x1)^3))];

% input force vector, B
B=[0 (10*sin(2*pi*50*t)) 0 0];

%output result as the equation 
Xdot=A+B';


end


