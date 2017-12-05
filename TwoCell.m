function [dy]=TwoCell(t,y,input,k1,m1,k2,m2)
% two unit cells, each with the double mass metamaterial system

% define X matrix, Xdot
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);
% first unit cell /\
%----------------------------------
% second unit cell \/
u1=y(5); % second mass 1 disp
u2=y(6); % second mass 1 velo
p1=y(7); % second mass 2 disp
p2=y(8); % second mass 2 velo
%----------------------------------
% input
f=input*(2*pi); %(Hz to rad/s)
%----------------------------------
% Define A matrix
A=[x2;...
    ((k2+2*k1)/m1)*x1-(k2/m1)*y1-(k1/m1)*u1;...
    y2;...
    -(k2/m2)*x1+(k2/m2)*y1;...
    u2;...
    ((k2+k1)/m1)*u1-(k2/m1)*p1-(k1/m1)*x1;...
    p2;...
    -(k2/m2)*u1+(k2/m2)*p1];

% Input excitation force
%sinusodial harmonic
H=0.1*(sin(f*t));

%sawtooth with random noise
% s=sawtooth(t);
% noise=awgn(s,0.5);   
% 
% % random noise
% R=0.1*randn(size(t));
B = [0 0 0 0 H/m1 0 0 0];


%output result as the equation 
dy=-A+B';
% dy=A*y;


end






