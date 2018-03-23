%% Plots the transient response of a forced 2DOF mass spring damper system 
clear all
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
%% Set the frequency range
f1=10;
f2=17;
f3=26;
f4=37;
f5=43;

freq_range=[f1 f2 f3 f4 f5]; % range from 10 Hz to 43 Hz 
% select which mass we want to observe the frequency response for 
res=1; % result of the m1 displacement = 1,
% m2 displacement = 3
freq_results=zeros(length(p:q),length(freq_range));% p and q are the new time start/end to cut 
% the steay state portion of the signal
% %% Total energy into the system 
% % need to find total energy into the system. We have frequencies, so need
% % to simulate for same length of time and calculate total energy under the
% % curve
% freq_input=zeros(length(p:q),length(freq_range));
% freq_input_int=zeros(5,1);
% 
% for i=1:length(freq_range)
%     input_signal=sin(2*pi*freq_range(i)*t(p:q));
%     % FFT loop
%     figure(i)
%     dt=abs(mean(diff(t)));  %average time step done
%     Fs=1/(dt);
%     % y = fft(ansys_amp_1);
%     % flop = (0:length(y)-1)*Fs/length(y);
%     n=length(t); %length of signal = number of samples
%     m=pow2(nextpow2(n));  %transform length
%     dft1=fft(input_signal,m); % DFT of signal
%     fr = (0:m-1)*(Fs/m);
%     fourier = abs(dft1);
%     f=Fs*(0:(n/2))/n;
%     freq1=fr(1:floor(m/2));
%     P1=fourier(1:floor(m/2));
%     plot(freq1,P1)
%     % plot(flop,abs(y),'LineWidth',2)
%     title(['FFT of AMM system at ', num2str(freq_range(i)), ' Hz'])
%     grid on
%     xlabel('Frequency,  (Hz)')
%     ylabel('|P1(f)|')
%     set(gca,'fontsize',20)
%     
%     % Integration loop
%     freq_input_int(i)=trapz(freq1,P1);
% end
%% Solve the model
for i=1:length(freq_range)
    omega=freq_range(i);
    options=odeset('InitialStep',dt,'MaxStep',dt);
    [t,result]=ode45(@(t,z) rhs(t,z,omega),t,z,options);
    
    % change result to show the steady state portion of the time history
    x=result(p:q,:); % x becomes the steady state result
    freq_results(:,i)=x(:,res); %store the displacement history in freq_results
end
toc
%% Results
% results are stored in freq_results vector
% integration loop of the different FFT results
freq_int=zeros(5,1);
freq_max=zeros(5,1);
% need to perform a for loop to create 5 figures of 5 FFTs.
for i=1:length(freq_range)
    
    % FFT loop
    figure(i)
    dt=abs(mean(diff(t)));  %average time step done
    Fs=1/(dt);
    % y = fft(ansys_amp_1);
    % flop = (0:length(y)-1)*Fs/length(y);
    n=length(t); %length of signal = number of samples
    m=pow2(nextpow2(n));  %transform length
    dft1=fft(freq_results(:,i),m); % DFT of signal
    fr = (0:m-1)*(Fs/m);
    fourier = abs(dft1);
    f=Fs*(0:(n/2))/n;
    freq1=fr(1:floor(m/2));
    P1=fourier(1:floor(m/2));
    plot(freq1,P1)
    % plot(flop,abs(y),'LineWidth',2)
    title(['FFT of AMM system at ', num2str(freq_range(i)), ' Hz'])
    grid on
    xlabel('Frequency,  (Hz)')
    ylabel('|P1(f)|')
    set(gca,'fontsize',20)
    
    % Integration loop
    freq_int(i)=trapz(freq1,P1);
    freq_max(i)=max(P1);
end

%% Plot int values as function of freq
figure
plot(freq_range,freq_int,'r--d')
title('Total energy within the system')
grid on
xlabel('Frequency,  (Hz)')
ylabel('Area under the curve')
set(gca,'fontsize',20)

%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=rhs(t,x,omega)
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