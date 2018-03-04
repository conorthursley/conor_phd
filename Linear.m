function [dy]=Linear(t,y,input,k1,m1,k2,m2)

% Define parameters of spring mass model
% k1=2e6;
% m1=1;
% k2=0.1*k1;
% m2=0.5*m1;


% define X matrix, Xdot
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);



% w0=sqrt(k2/m2);
func=input(1)*sin(input(2)*t);

% Define A matrix
A=[x2; ...
    1+x1*((k1+k2)/m1)-y1*((k2/m1)); ...
    y2; ...
    -((k2/m2)*x1)+((k2/m2)*y1)];

% Input excitation force
%sinusodial harmonic
H=input(2)*sin(input(1));

%sawtooth with random noise
% s=sawtooth(t);
% noise=awgn(s,0.5);   
% 
% % random noise
% R=input(1)*randn(size(t));
% B = [0 1 0 0];


%output result as the equation 
dy=A;
% dy=A*y;

end






