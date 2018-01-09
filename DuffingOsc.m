function [dy]=DuffingOsc(t,y,input,k1,m1,k2,m2)

%-----------------------------------
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);
%-----------------------------------
f=input(1)*(2*pi); 
k3=10*k2;

A=[x2;...
    -((k1/m1)+(k2/m1))*(x1)+((k3/m1)*y1^3); ...
    y2;...
    -((k3/m2)*(y1-x1))-((k3/m2)*((y1-x1)^3))];
    
% Input excitation force
%sinusodial harmonic
H=input(1)*(sin(f*t));
B = [0 H/m1 0 0];


%output result as the equation 
dy=A-B';
% dy=A*y;

end
