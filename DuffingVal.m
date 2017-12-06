function [dy]=DuffingVal(~,y,~,k1,m1,k2,~)

%-----------------------------------
x1=y(1);
x2=y(2);
%-----------------------------------
% f=input(1)*(2*pi); 
d=0.01; %damping
A=[x2;-(d/m1)*x2-((k1/m1)*x1+(k2/m1))*(x1.^3)];
    
% A=y(2);
% B=gamma*cos(omega*2*pi*t)-delta*y(2)-beta*y(1)-alpha*y(1).^3;
% dy=[A; B];


% Input excitation force
%sinusodial harmonic
% H=(input(2)*(sin(f*t)));

B=[0 1];

%output result as the equation 
% dy=A*y;
dy=A+B';
end

