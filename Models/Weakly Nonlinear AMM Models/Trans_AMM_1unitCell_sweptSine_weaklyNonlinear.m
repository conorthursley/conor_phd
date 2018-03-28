%% Plots the transient response of a forced 2DOF mass spring damper system 
clear 
tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;    % [s] delta t
% for loop parameters
t_end=800;   % t limit
t=0:dt:t_end;      % [s] time scale
t_find=400; % the time to safely assume SS has been reached 600 seconds after initial transient begins
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
freq_step=0.1;
swept_sine_range=10:freq_step:150; % range from 10 Hz to 43 Hz in steps of 0.25 Hz

amplitude=zeros(length(swept_sine_range),4); %vector to store amps of displacement and velocity 

%% set the nonlinear strength
sigma=[5 10 25 50 100 250 500 750];

%% Solve the model
% for j=1:length(sigma)
    k3=sigma(6);
    for i=1:length(swept_sine_range)
        omega=swept_sine_range(i);
        options=odeset('InitialStep',dt,'MaxStep',dt);
        [t,result]=ode45(@(t,z) rhs(t,z,omega,k3),t,z,options);
        
        % change result to show the steady state portion of the time history
        x=result(p:q,:); % x becomes the steady state result
        amplitude(i,:)=max(x);
    end
toc
%% Results
% extract displacement amplitudes from vector, "amplitude"
% amplitude = [displacement1 velo1 disp2 velo2]
% displacement
m1_disp=amplitude(:,1);
m2_disp=amplitude(:,3);
% velocity 
m1_velo=amplitude(:,2);
m2_velo=amplitude(:,4);

%% Plot the results
% Plot the time series
figure
plot1=loglog(swept_sine_range,m1_disp,swept_sine_range,m2_disp);
set(plot1,'LineWidth',2)
xlabel('Freqency range (Hz)'); ylabel('Displacment (m)');
title(['Swept Sine results for displacement with k3 = ',num2str(sigma(j))])
grid on
legend 'mass1' 'mass2' 
set(gca,'fontsize',20) 

% Plot the Kinetic Energy ratio
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
title(['Swept Sine results for energy distribution rate of the AMM model with k3 = ',num2str(sigma(j))])
grid on
legend_text=['\theta=',num2str(theta)];
legend(legend_text,'Analytical Result','FontAngle','italic','Interpreter','Latex')
set(gca,'fontsize',20) 
% end

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