%% setup the freqeuncy range 17 to 22 in steps of 0.25
clear all
lt=300; %length of time
il=2.5; %interval length

low=14; %low interval
high=24; %high interval

hf=high; %highest frequency
dt=1/(10*hf); %time step size
Fsteps=(high-low)/(lt/il);
freq=low+Fsteps:Fsteps:high;
step_freq=low;


%% plotting duffing relation
% 
%parameters
m=0.1;
a=1000/m;
y=227/m;
d=0;
B=10/m;
n=1;
dis=[];
 %linear oscillator - as B is the nonlinear term
% z is steady state displacment and w is frequency
for ii=1:(lt/il) %frequency intervals
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
%%
file='C:\ANSYS\Temp\Validation\DuffingValDec17\DuffingNumerical.csv';
csvwrite(file,dis);

%% figure
figure(1)
plot(freq,dis(:,1),'b',freq,dis(:,2),'g',freq,dis(:,6),'k')
grid on
xlabel 'Input Frequency, Hz'
ylabel 'Displacement, mm'
title 'Displacement of U1 and U2 as a function of w'
legend S1 S2 S3
% figure
% xlabel('frequency, w','FontSize',14) % x-axis label
% ylabel('amplitude, z','FontSize',14) % y-axis label
% title('Duffing Relation $\ddot{x}$ + $\delta\dot{x}$ + $\alpha{x}$ $\beta{x}^3$ = $\gamma{cos(\omega{t})}$', 'Interpreter','latex' ,'FontSize',14)
% legend({'\beta = %d','\beta = %d','\beta = %d','\beta = %d','\beta = %d'},'FontSize',14)
% grid on


%% Old Script 
%% plotting duffing relation
% 
% close all
m=0.1;
a=25.6/m;
y=200/m;
d=0;
B=0.256/m;
%parameters
b=[0 B*10 B B/10]; %linear oscillator - as B is the nonlinear term
% z is steady state displacment and w is frequency
% figure
figure(2)
for i=1:length(b)
   B=b(i);
   
   f = @(w,z) ((w.^2-a-(3/4)*B.*z.^2)^2+(d.*w).^2).*z.^2-y^2;

   fimplicit(f,[0 150 0 50])
   hold on
end

xlabel('frequency, w','FontSize',14) % x-axis label
ylabel('amplitude, z','FontSize',14) % y-axis label
title('Duffing Relation $\ddot{x}$ + $\delta\dot{x}$ + $\alpha{x}$ $\beta{x}^3$ = $\gamma{cos(\omega{t})}$', 'Interpreter','latex' ,'FontSize',14)
legend({'\beta = %d','\beta = %d','\beta = %d','\beta = %d','\beta = %d'},'FontSize',14)
grid on
