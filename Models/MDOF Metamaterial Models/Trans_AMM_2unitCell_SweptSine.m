%% Plots the transient response of a forced 2 unit cell MDOF metamaterial model
clear all
tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;    % [s] delta t
% for loop parameters
t_end=1000;   % t limit
t=0:dt:t_end;      % [s] time scale
t_find=600; % the time to safely assume SS has been reached 600 seconds after initial transient begins
p=find(t==600); q=find(t==t_end);

mass1=101.10e-3;		% [kg]
mass2=46.47e-3;
stiff1=117;    % [N/m]
stiff2=37*2;

w2=sqrt(stiff2/mass2)/(2*pi);
theta=mass2/mass1;


%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
z=zeros(1,2*4); % n=2 and there are 4 DOF per unit cell, hence 4 initial conditions;
%% Set the frequency range
freq_step=0.1;
swept_sine_range=0:freq_step:13; % range from 10 Hz to 45 Hz in steps of 0.25 Hz

amplitude=zeros(length(swept_sine_range),length(z)); %vector to store amps of displacement and velocity 

%% Solve the model
parfor i=1:length(swept_sine_range)
    t=0:dt:t_end;
    omega=swept_sine_range(i);
    options=odeset('InitialStep',dt,'MaxStep',dt);
    [t,result]=ode45(@(t,z) rhs(t,z,omega),t,z,options);
    
    % change result to show the steady state portion of the time history
    x=result(p:q,:); % x becomes the steady state result
    amplitude(i,:)=max(x);
end
t_new=t(p:q);
toc
%% Results
% extract displacement amplitudes from vector, "amplitude"
% amplitude = [displacement1 velo1 disp2 velo2 disp'n' velo'n' .....]
%----------First unit cell--------------
% displacement
m1_disp=amplitude(:,1);
m2_disp=amplitude(:,3);
% velocity 
m1_velo=amplitude(:,2);
m2_velo=amplitude(:,4);
%----------Last unit cell---------------
% displacement, 'last mass' - 'lm'
lm1_disp=amplitude(:,end-3);
lm2_disp=amplitude(:,end-1);
% velocity 
lm1_velo=amplitude(:,end-2);
lm2_velo=amplitude(:,end);

%% Plot the results
% Plot the time series
figure
plot1=plot(swept_sine_range,mag2db(m1_disp),swept_sine_range,mag2db(m2_disp)...
    ,swept_sine_range,mag2db(lm1_disp),swept_sine_range,mag2db(lm2_disp),data(:,2),data(:,1),'r',data1(:,2),data1(:,1),'g');
set(plot1,'LineWidth',2)
xlabel('Freqency range (Hz)'); ylabel('Transmittance (dB)');
title('Swept Sine results for displacement')
grid on
legend 'mass1' 'mass2' 'last mass1' 'last mass2'
set(gca,'fontsize',20) 
%% Transmission
% U=Xn/X1
U=abs(lm2_disp)./(m2_disp);
figure
plot(swept_sine_range,20*log10((U)),'b',data(:,2),data(:,1),'r',data1(:,2),data1(:,1),'g')
xlabel('Freqency range (Hz)'); ylabel('Transmittance (dB)');
title('Transmittance of the finite system')
grid on
set(gca,'fontsize',20) 
% periodogram(x(:,3))
%% Plot the Kinetic Energy ratio
%------------------numerical
KE1=0.5*mass1.*(m1_velo.^2);
KE2=0.5*mass2.*(m2_velo.^2);
RDR=KE2./(KE1+KE2);   %ratio
%------------------analytical
RDR_a=theta./((1-(swept_sine_range./w2).^2).^2+theta);

figure
plot1=plot(swept_sine_range/w2,RDR); hold on
plot2=plot(swept_sine_range/w2,RDR_a,'r--');
set([plot1 plot2],'LineWidth',3.5)
xlabel('Normalised freqency \omega/\omega_2'); ylabel('Ratio of Kinetic Energy');
title('Swept Sine results for energy distribution rate of the AMM model')
grid on
legend_text=['\theta=',num2str(theta)];
legend(legend_text,'Analytical Result','FontAngle','italic','Interpreter','Latex')
set(gca,'fontsize',20) 
%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=rhs(t,x,omega)
        mass1=101.10e-3;		% [kg]
        mass2=46.47e-3;
        stiff1=117;    % [N/m]
        stiff2=37*2;
        damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.002;
        f=1; %*(stepfun(t,0)-stepfun(t,0.01));
        w=omega; % Hz, forcing frequency 
     
        %---------------------------------------
        % first unit cell
        % first mass
        dxdt_1 = x(2);
        dxdt_2 = -((2*damp1+damp2)/mass1)*x(2)- ((2*stiff1+stiff2)/mass1)*x(1) + (stiff1/mass1)*x(5)+ (damp1/mass1)*x(6)...
            +(stiff2/mass1)*x(3)+(damp2/mass1)*x(4)+((f/mass1)*sin(2*pi*w*t)*stiff1);
        % second mass
        dydt_1= x(4);
        dydt_2= -(stiff2/mass2)*x(3) - (damp2/mass2)*x(4) + (stiff2/mass2)*x(1) + (damp2/mass2)*x(2);
        %---------------------------------------
        % Last unit cell (cell 2)
        % first mass
        dxdt_3 = x(6);
        dxdt_4 = -((2*damp1+damp2)/mass1)*x(6)- ((2*stiff1+stiff2)/mass1)*x(5)...
            + (stiff1/mass1)*x(1)+ (damp1/mass1)*x(2) +(stiff2/mass1)*x(7)+(damp2/mass1)*x(8);
        % second mass
        dydt_3= x(8);
        dydt_4= -(stiff2/mass2)*x(7) - (damp2/mass2)*x(8) + (stiff2/mass2)*x(5) + (damp2/mass2)*x(6);
        %---------------------------------------
        %% -----------Final Solution-------------
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2;dxdt_3;dxdt_4;dydt_3;dydt_4];
end