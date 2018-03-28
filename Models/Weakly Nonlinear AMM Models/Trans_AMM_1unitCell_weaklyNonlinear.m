%% Plots the transient response of a forced 2DOF mass spring damper system 
clear all
tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;    % [s] delta t
% for loop parameters
t_end=2000;   % t limit
t=0:dt:t_end;      % [s] time scale
t_find=600; % the time to safely assume SS has been reached 600 seconds after initial transient begins
p=find(t==600); q=find(t==t_end);
% model parameters
mass1=0.1;		% [kg]
mass2=mass1*0.5;
stiff1=1000;    % [N/m]
stiff2=1.5*stiff1;
% nonlinear parameter
stiff3=0*stiff2;
k3=stiff3;
% driving frequency 
omega=10;


w2=sqrt(stiff2/mass2)/(2*pi);
theta=mass2/mass1;

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0e-3;
initial_dxdt = 0;
initial_y    = 0e-3;
initial_dydt = 0;

z=[initial_x initial_dxdt initial_y initial_dydt];

%% Solve the model

options=odeset('InitialStep',dt,'MaxStep',dt);
[t,result]=ode45(@(t,z) rhs(t,z,omega,k3),t,z,options);
x=result(p:q,:); % x becomes the steady state result
t_new=t(p:q);
toc
%% Plot the results
% Plot the time series
% Original time series
figure
plot1=plot(t,result(:,1),t,result(:,3));
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' 
set(gca,'fontsize',20) 
% Clipped time series 
figure
plot1=plot(t_new,x(:,1),t_new,x(:,3));
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' 
set(gca,'fontsize',20) 
%% Plot the phase plane
% x=result;
figure
plot1=plot(x(:,2),x(:,4));
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Phase Portrait')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' 
set(gca,'fontsize',20) 
%% FFT
figure
dt=abs(mean(diff(t)));  %average time step done
Fs=1/(dt);
% y = fft(ansys_amp_1);
% flop = (0:length(y)-1)*Fs/length(y);
n=length(t); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(x(:,1),m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft1);
f=Fs*(0:(n/2))/n;
freq1=fr(1:floor(m/2));
P1=fourier(1:floor(m/2));
plot(freq1,P1)
% plot(flop,abs(y),'LineWidth',2)
title(['FFT of AMM system at ', num2str(omega), ' Hz'])
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',20)

%% Work and Energy functions
%-----------Results----------------------
% extract displacement amplitudes from vector, "x"
% x = [displacement1 velo1 disp2 velo2]
% displacement
m1_disp=x(:,1);
m2_disp=x(:,3);
% velocity 
m1_velo=x(:,2);
m2_velo=x(:,4);
%------------Work-----------------------------
% Looking at the dynamic model, apply the energy approach for KE and PE and
% find the total KE and PE of the system
% KE=KE1 + KE2
% PE=PEu1+PEu1+PE(u2-u1) (PEu1 is done twice as we have two k1 springs on
% either side of the model)
% KE------------
% Kinetic Energy = 0.5*m_i*v_i^2
KE=0.5*mass1*m1_velo.^2 + 0.5*mass2*m2_velo.^2;
% PE------------
% Potential Energy = 0.5*m_i*v_i^2
PE=stiff1*m1_disp.^2+0.5*stiff2*(m2_disp-m1_disp).^2;

figure
plot(t_new,KE,'r',t_new,PE,'b',t_new,(PE+KE),'k:');
xlabel('t'); ylabel('Work/PE/KE');
title(['Work/Energy Calculations at ', num2str(omega), ' Hz'])
grid on
legend 'KE' 'PE' 'Work' 
set(gca,'fontsize',20) 

%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=rhs(t,x,omega,k3)
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