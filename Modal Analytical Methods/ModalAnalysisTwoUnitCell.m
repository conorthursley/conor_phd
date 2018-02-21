%% Modal analysis 
% steps
% 1) solve undamped eigenvalue problem for resonant frequencies and mode
% shapes (eigenvalues and eigenvectors) - useful for understanding basic
% motion of the system

%---------------------------
% Parameters
%---------------------------

k1=1000;
k2=1.5*k1;
m1=0.1;
m2=0.5*m1;
%---------------------------
cr1=2*sqrt(k1*m1);
cr2=2*sqrt(k2*m2);
zeta=0.0002;
c1=cr1*zeta;
c2=cr2*zeta;

Ug=0.1; % Forcing amplitude
%---------------------------
% Matrices
%---------------------------
n=2; % no. of cells
m_vector=[m1 m2];
m_array=repmat(m_vector,n);
M=diag(m_array(1,:));
%---------------------------
K=[2*k1+k2 -k2 -k1 0;-k2 k2 0 0;-k1 0 2*k1+k2 -k2;0 0 -k2 k2];

C=[2*c1+c2 -c2 -c1 0; -c2 c2 0 0;-c1 0 2*c1+c2 -c2; 0 0 -c2 c2];

%---------------------------
% Need to create system matrix A to solve for eigvalues
%---------------------------

% by inspection;
A=[0 1 0 0 0 0 0 0; -(2*k1+k2)/m1 2*c1+c2/m1 k2/m1 c2/m1 k1/m1 c1/m1 0 0; 0 0 0 1 0 0 0 0;k2/m2 c2/m2 -k2/m2 -c2/m2 0 0 0 0;...
    0 0 0 0 0 1 0 0;k1/m1 c1/m1 0 0 -(2*k1+k2)/m1 2*c1+c2/m1 k2/m1 c2/m1; 0 0 0 0 0 0 0 1;0 0 0 0 k2/m2 c2/m2 -k2/m2 -c2/m2];

B=[0; k1*Ug; 0;0;0;0;0;0];
CC=[1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0;0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0];
D=[0];
sys=ss(A,B,CC,D);
bode(sys)
% eigenvalues
[Zm,lambda]=eig(A);

lambdad=diag(lambda);
[lambdao,index]=sort(abs(imag(lambdad)));

% reorder eig values from low to high

w1=lambdao(1);
w2=lambdao(3);
w3=lambdao(5);
w4=lambdao(7);
%---------------------------
% Modal Matrix
%---------------------------
% will need to do this for different cases, so write code to determine it
% pls.

% determine mode shapes of each natural frequency
a1=-m2/k2*w1^2+1;
a2=-m2/k2*w2^2+1;

phi_normal=[a1 a2; 1 1];
% 2) use eigenvectors to uncouple or diagonalies the original set of
% coupled equations, allowing the solution of uncoupled SDOF problems
% instead of solving a set of coupled equations
%---------------------------
% Mass Normalised Modal Matrix
%---------------------------
den_mass1 = sqrt(m1*a1^2+m2);
den_mass2 = sqrt(m1*a2^2+m2);

phi=[a1/den_mass1 a2/den_mass2; 1/den_mass1 1/den_mass2];
%---------------------------
% Mass Normalised Mass Matrix

Mn=phi'*M*phi;

% Mass Normalised Stiffness Matrix

Kn=phi'*K*phi;

% Mass Normalised damping Matrix
Cn=phi'*C*phi;
% Mass Normalised Force vector

f=[k1*Ug;0]; %force input - usually multiplie by sin(wt)
P=phi'*f;
%---------------------------

%% uncoupled steady state equations of motion
% input freuency
w_bar=50; %rad.s
%---------------------------
eta1=w_bar/w1;
eta2=w_bar/w2;
gamma1=atan((-2*zeta*eta1)/(1-eta1^2));
gamma2=atan((-2*zeta*eta2)/(1-eta2^2));
den1=(sqrt((1-eta1^2)^2)+(2*zeta*eta1)^2);
den2=(sqrt((1-eta2^2)^2)+(2*zeta*eta2)^2);
%---------------------------
num1=(P(1)/w1^2);
num2=(P(2)/w2^2);
Z1=num1/den1;
Z2=num2/den2;

t=linspace(0,1,1000);
lem=sin(w_bar*t);

%---------------------------
% plots
u1=(Z1/w1)*lem.*sin(w1*t);
u2=(Z2/w2)*lem.*sin(w2*t);

figure
plot(t,u1,'b',t,u2,'r')

% tfestimate(lem,u2)
