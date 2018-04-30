%% Comparison of Linear and nonlinear systems 1U system 
% Plots the transient response of a forced 2DOF mass spring damper system 
clear 
tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;    % [s] delta t
% for loop parameters
t_end=1500;   % t limit
t=0:dt:t_end;      % [s] time scale
t_find=600; % the time to safely assume SS has been reached 600 seconds after initial transient begins
p=find(t==600); q=find(t==t_end);

mass1=0.1;		% [kg]
mass2=mass1*0.5;
stiff1=1000;    % [N/m]
stiff2=1.5*stiff1;

w2=sqrt(stiff2/mass2)/(2*pi);
theta=mass2/mass1;

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0e-3;
initial_dxdt = 0;
initial_y    = 0e-3;
initial_dydt = 0;

z=[initial_x initial_dxdt initial_y initial_dydt];

stiff3=stiff2*1600;
k3=stiff3;

omega=20;
%% Solve both of the models

t=0:dt:t_end;      % [s] time scale
options=odeset('InitialStep',dt,'MaxStep',dt);
%linear case
[t1,result1]=ode45(@(t,z) Lrhs(t,z,omega),t,z,options);
% change result to show the steady state portion of the time history
x1=result1(p:q,[1,3]); % x becomes the steady state result
%nonlinear case
[t2,result2]=ode45(@(t,z) NLrhs(t,z,omega,k3),t,z,options);
x2=result2(p:q,[1,3]); % x becomes the steady state result
toc
%% Results
% new time
t1_new=t1(p:q);
t2_new=t2(p:q);

%% Plot disp values
figure
plot(t1_new,x1,'g',t2_new,x2,'r')
title('Displacement')
grid on
xlabel('Time,  (Hz)')
ylabel('Amplitude')
% legend 'linear' 'nonlinear'
set(gca,'fontsize',20)
%% FFT
figure
%linear
subplot(2,1,1)
dt=abs(mean(diff(t1_new)));  %average time step done
Fs=1/(dt);
n=length(t1_new); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(x1(:,1),m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft1);
f=Fs*(0:(n/2))/n;
freq1=fr(1:floor(m/2));
P1=fourier(1:floor(m/2));
plot(freq1,P1)
title(['FFT of AMM system at ', num2str(omega), ' Hz'])
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',20)
%Nonlinear
subplot(2,1,2)
dt=abs(mean(diff(t2_new)));  %average time step done
Fs=1/(dt);
n=length(t2_new); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(x2(:,1),m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft1);
f=Fs*(0:(n/2))/n;
freq1=fr(1:floor(m/2));
P1=fourier(1:floor(m/2));
plot(freq1,P1)
title(['FFT of AMM system at ', num2str(omega), ' Hz'])
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',20)
%% Nonlinear Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=NLrhs(t,x,omega,k3)
        mass1=0.1;		% [kg]
        mass2=mass1*0.5;
        stiff1=1000;    % [N/m]
        stiff2=1.5*stiff1;
        stiff3=k3;
        damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.002;
        f=1; %*(stepfun(t,0)-stepfun(t,0.01));
        w=omega; % Hz, forcing frequency 
        u=x(1);    %disp mass2
        du=x(2);    %velo mass1
        v=x(3);   %disp mass2
        dv=x(4);  % velo mass2
     
        %---------------------------------------
        % first unit cell
        % first mass
        dxdt_1 = du;
        dxdt_2 = -((2*damp1+damp2)/mass1)*du- ((2*stiff1)/mass1)*u-(stiff2/mass1)*(u-v) -...
            (stiff3/mass1)*(u-v)^3+(damp2/mass1)*dv+(f/mass1)*sin(2*pi*w*t);
        % second mass
        dydt_1= dv;
        dydt_2= -(stiff2/mass2)*(v-u)-(stiff3/mass2)*(v-u)^3 - (damp2/mass2)*dv + (damp2/mass2)*du;
        %---------------------------------------
                
        % final solution 
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2];
end
%% Linear Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=Lrhs(t,x,omega)
        mass1=0.1;		% [kg]
        mass2=mass1*0.5;
        stiff1=1000;    % [N/m]
        stiff2=1.5*stiff1;
        damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.002;
        f=1; %*(stepfun(t,0)-stepfun(t,0.01));
        w=omega; % Hz, forcing frequency 
     
        %---------------------------------------
        % first unit cell
        % first mass
        dxdt_1 = x(2);
        dxdt_2 = -((2*damp1+damp2)/mass1)*x(2)- ((2*stiff1+stiff2)/mass1)*x(1) +(stiff2/mass1)*x(3)+(damp2/mass1)*x(4)...
          +(f/mass1)*sin(2*pi*w*t);
        % second mass
        dydt_1= x(4);
        dydt_2= -(stiff2/mass2)*x(3) - (damp2/mass2)*x(4) + (stiff2/mass2)*x(1) + (damp2/mass2)*x(2);
        %---------------------------------------
                
        % final solution 
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2];
end
