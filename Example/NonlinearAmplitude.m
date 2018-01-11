%% Validation of nonlinear mass spring system of Duffing type
% can also try van der pol oscillator later on
%% File read from APDL simulation


file = 'U:\_PhD\APDL\Validation\DuffingValDec17\DuffOneUnitTrans174.csv';
M=csvread(file,1,0); %start reading from row 3, column 1

time = M((1:end),1);      % excitation frequency
amp =abs((M((1:end),4)));
w0=(zeros(length(time),1))';
w1=(zeros(length(time),1))';
%-------------------------------------------------------
% setup the freqeuncy range 17 to 22 in steps of 0.25
lt=300; %length of time
il=10; %interval length
A=2000; %amplitude 

low=1.5; %low interval
high=60; %high interval

hf=high; %highest frequency
dt=1/(10*hf); %time step size
Fsteps=(high-low)/(lt/il);
freq=low+Fsteps:Fsteps:high;
ranger=1:(length(M)/(lt/il));
step_freq=low;
for ii=1:(lt/il)-1 %frequency intervals
    w0(:,ranger)=step_freq;
    mu=mean(amp(ranger));
    w1(:,ranger)=mu;   % average 
    step_freq=step_freq+Fsteps;
    ranger=ranger+round((length(M)/(lt/il)));
end
figure
hold on
plot(w0,(w1),'c-*')
grid on
xlabel ('Input Frequency, Hz','FontSize',18)
ylabel ('Displacement, mm','FontSize',18)
title ('Displacement of Duffing Oscillator with increasing and decreasing frequency','FontSize',18)
legend S1 S2 S3 var1 var2 var3 var4


%% plotting duffing relation
% 
%parameters
m=0.5;
a=1000/m;
y=200/m;
d=0;
B=10/m;
n=1;
dis=[];
%[-0.5 -0.05 0 0.05 0.5]; %linear oscillator - as B is the nonlinear term
% z is steady state displacment and w is frequency
for J=1:(lt/il) %frequency intervals
   w=step_freq;
    w=w*2*pi;
    syms z
    equ = ((w.^2-a-(3/4)*B*z^2)^2+(d*w).^2)*z^2-y^2;
    var1=solve(equ,z);
    S=double(var1);
    dis(n,:)=[S];
    
    % counters
    n=n+1;
    step_freq=step_freq+Fsteps;
end


%% Plots and plotting
 % Plot of U1 and U2 displacements against input frequency
% figure
hold on
plot(w0,(w1),'r-o',freq,dis(:,1),'b',freq,dis(:,2),'g',freq,dis(:,6),'k') %,freq,dis(:,3),'k')
grid on
xlabel 'Input Frequency, Hz'
ylabel 'Displacement, mm'
title 'Displacement of U1 and U2 as a function of w'
