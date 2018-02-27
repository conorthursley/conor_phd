%% First Order Analytical approximate displacements
%----------------------
%Analysis of nonlinear vibration of coupled systems with cubics
%----------------------
% variables
k1=1e5;
<<<<<<< HEAD
k2=1e4;
=======
k2=1e7;
>>>>>>> master
m=1;

alpha=k1/m; %for a single mass value for both masses
beta=k2/m;
A=1; %amplitude, I think. leaving as 1 for now
t=0:0.005:1; %time range
f=2; %Hz
w=linspace(0,2.2,length(t)); %forcing frequency
zeta=0.05;
B=3; 

M=1./sqrt((w.^2-1-(3/4)*B*1).^2+(2*zeta*w).^2);
plot(w,M)
grid


%% Displacement of mass 1, u(t)

% u=-(1./(9*w))*A.*cos(w.*t).*(9*alpha+6*beta*A^2+A*beta*(cos(w.*t).^2));
% v=u+A.*cos(w.*t);
% plot(t,u,t,v)