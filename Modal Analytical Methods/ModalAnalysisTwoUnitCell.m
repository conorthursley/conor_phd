%% Modal analysis Two Unit Cell
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
K=[2*k1+k2 -k2 -k1 0;-k2 k2 0 0;-k1 0 k1+k2 -k2;0 0 -k2 k2];

C=[2*c1+c2 -c2 -c1 0; -c2 c2 0 0;-c1 0 2*c1+c2 -c2; 0 0 -c2 c2];

%---------------------------
% Need to create system matrix A to solve for eigvalues
%---------------------------

% by inspection; and with no damping
A1=[0 1 0 0 0 0 0 0; -(2*k1+k2)/m1 0 k2/m1 0 k1/m1 0 0 0; 0 0 0 1 0 0 0 0;k2/m2 0 -k2/m2 0 0 0 0 0;...
    0 0 0 0 0 1 0 0;k1/m1 0 0 0 -(k1+k2)/m1 0 k2/m1 0; 0 0 0 0 0 0 0 1;0 0 0 0 k2/m2 0 -k2/m2 0];
% by inspection; and with damping
A2=[0 1 0 0 0 0 0 0; -(2*k1+k2)/m1 2*c1+c2/m1 k2/m1 c2/m1 k1/m1 c1/m1 0 0; 0 0 0 1 0 0 0 0;k2/m2 c2/m2 -k2/m2 -c2/m2 0 0 0 0;...
    0 0 0 0 0 1 0 0;k1/m1 c1/m1 0 0 -(2*k1+k2)/m1 2*c1+c2/m1 k2/m1 c2/m1; 0 0 0 0 0 0 0 1;0 0 0 0 k2/m2 c2/m2 -k2/m2 -c2/m2];
B=[0; k1*Ug; 0;0;0;0;0;0];
CC=[1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0;0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0];
D=[0];
sys=ss(A1,B,CC,D);
bode(sys)
% eigenvalues
[Xm,lambda]=eig(A1);

lambdad=diag(lambda);
[lambdao,index]=sort(abs(imag(lambdad)));

% reorder eig values from low to high
W=zeros(4,1);
W(1)=lambdao(1);
W(2)=lambdao(3);
W(3)=lambdao(5);
W(4)=lambdao(7);

%% ---------------------------
% Modal Matrix
%---------------------------
% will need to do this for different cases, so write code to determine it
% pls.
% write symbolic matrix of A and Z vector to solve for the eigenvector
% ratios to get the modal matrix
% Zm=[z11 z12 z13 z14|
%    |z21 z22 z23 z24|
%    |z31 z32 z33 z34|
%    |z41 z42 z43 z24]
%     Z1  Z2  Z3  Z4    Vectors
% for z2 is chosen as a reference
% thus z1/z2 = (k2-w1^2*m2)/k2
%      z3/z2 = (((2*k1+k2-w^2*m1)*(k2-w^2*m2))/(k2*k1))-(k2/k1)
%      z4/z2 = ((2*k1+k2-w^2*m1)*(z3/z2)-k1)/k2)
% for z1 as a reference (fixed-fixed)
% thus z2/z1 = k2/(k2-w1^2*m2)
%      z3/z1 = ((2*k1+k2-w^2*m1)-k2^2/(k2-w1^2*m2))/k1
%      z4/z1 = ((2*k1+k2-w^2*m1)/k2)*(z3/z1)-(k1/k2)
% for z1 as a reference (free-free)
% thus z2/z1 = k2/(k2-w1^2*m2) = an
%      z3/z1 = ((k1+k2-w^2*m1)-k2*(an))/k1 = bn
%      z4/z1 = (((k1+k2-w^2*m1)*bn)-k1)/k2

Z1=ones(1,4);
Z2=zeros(1,4);
Z3=zeros(1,4);
Z4=zeros(1,4);
for i=1:4
    Z2(i)=k2/(k2-W(i)^2*m2);
    Z3(i)=((2*k1+k2-W(i)^2*m1)-k2*(Z1(i)))/k1;
    Z4(i)=(((k1+k2-W(i)^2*m1)*Z3(i))-k1)/k2;
end
Zm=[Z1;Z2;Z3;Z4];
Q=diag(M);
q1=sqrt((Q(1)*(Zm(1,1))^2+Q(2)*(Zm(2,1))^2+Q(3)*(Zm(3,1))^2+Q(4)*(Zm(4,1))^2));
q2=sqrt((Q(1)*(Zm(1,2))^2+Q(2)*(Zm(2,2))^2+Q(3)*(Zm(3,2))^2+Q(4)*(Zm(4,2))^2));
q3=sqrt((Q(1)*(Zm(1,3))^2+Q(2)*(Zm(2,3))^2+Q(3)*(Zm(3,3))^2+Q(4)*(Zm(4,3))^2));
q4=sqrt((Q(1)*(Zm(1,4))^2+Q(2)*(Zm(2,4))^2+Q(3)*(Zm(3,4))^2+Q(4)*(Zm(4,4))^2));

Zn1=Zm(:,1)/q1;
Zn2=Zm(:,2)/q2;
Zn3=Zm(:,3)/q3;
Zn4=Zm(:,4)/q4;

Zn=[Zn1 Zn2 Zn3 Zn4];

Zn'*M*Zn
Zn'*K*Zn
%% 
z21=(k2-w1^2*m2)/k2;
z31=(w1^2*(m2-m1)+2*k1)/k1;
z44=((2*k1+k2-w4^2*m1)*(w4^2*(m2-m1)+2*k1)-k1^2)/(k1*k2);
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
