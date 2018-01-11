function Xdot=nonLinear(t,y,input,NLspring,k1,m1,k2,m2)
%Define parameters from model
% k1=1000; % Not Used
%------------------------------------------
%% check this with the larger values of m2 to show Carl and Luke

%-----------------------------------------------------------
%% 
% m1=1;
% k2=100;
% m2=500*m1;
% k3=2e4;
% D=1e10; %damping coefficient 

%define 'a' coefficient to resolve large solution issue
% a=1e-6;
%define matrices from initial conditions
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);
%----------------------------------
%set nonlinear term, k3, to 0
% NL=0;%********nonlinear ratio******** nonlinear/linear
% k3=k2*NL; 
% H=input(1)*cos(2*pi*input(2)*t); 
H=input(t);
% value = [double(y(3)<=(y(1)-(5e-4))); double(y(3)>=y(1)+(5e-4))]
%define A matrix for EOM
A=[x2; ...
    H/m1-((k1/m1))*x1-(NLspring(x1-y1))/m1
%     ((k2/m1))*(x1-y1)-((k3/m1)*((x1-y1).^3)); ...
    y2;...
    -(NLspring(y1-x1))/m1];
% ((k2/m2)*(y1-x1))-((k3/m2)*((y1-x1).^3))
%----------------------------------
% input excitation
% H=((input(1))*(sin(2*pi*input(2)*t)));
% white noise
% H=input(1)*sin(2*pi*input(2)*t); %*randn(1,t);
% 
% % input force vector, B
% B = [0 H/m1 0 0];
%----------------------------------
%output result as the equation 
Xdot=A;


end



