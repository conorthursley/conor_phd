%% solving the duffing equation theoretically to get our amplitude vs frequency curves

clear all
% close all
m=1;
kl=1e+3; %linear stiffness, k2
kn=100; %nonlinear stiffness
c=1; %damping
f=200; %forcing amplitude 
%-------------------------------------------------
w0=sqrt(kl/m);
% omega=w/w0;
syms u
kms=[];
disp=[];
n=1;
for w=4:0.1:8
    kms(n,1)=w;
    equ=31.6228+1.1859*(u^2)+sqrt((-1+10/((u^2))))-w;  
    S=solve(equ,u); %,'MaxDegree',3);
    dis=double(S);
    disp(n,1)=dis(2)
    n=n+1;
end

%-------------------------------------------------