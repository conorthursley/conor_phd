%% Plots the transient response of a forced 10 unit cell MDOF metamaterial model
clear 
tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;    % [s] delta t
% for loop parameters
t_end=1000;   % t limit
t=0:dt:t_end;      % [s] time scale
t_find=600; % the time to safely assume SS has been reached 600 seconds after initial transient begins
p=find(t==600); q=find(t==t_end);

mass1=0.1;		% [kg]
mass2=0.5*mass1;
stiff1=1000;    % [N/m]
stiff2=1500;

w2=sqrt(stiff2/mass2)/(2*pi);
theta=mass2/mass1;


%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
z=zeros(1,5*4); % n=5 and there are 4 DOF per unit cell, hence 4 initial conditions;
%% Set the frequency range
freq_step=0.25;
swept_sine_range=freq_step:freq_step:60; % range from 10 Hz to 45 Hz in steps of 0.25 Hz

amplitude=zeros(length(swept_sine_range),length(z)); %vector to store amps of displacement and velocity 

%% set the nonlinear strength
k3=1600*stiff2;
%% Solve the model
parfor i=1:length(swept_sine_range)
    t=0:dt:t_end;
    omega=swept_sine_range(i);
    options=odeset('InitialStep',dt,'MaxStep',dt);
    [t,result]=ode45(@(t,z) rhs(t,z,omega,k3),t,z,options);
    
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
    ,swept_sine_range,mag2db(lm1_disp),swept_sine_range,mag2db(lm2_disp));
set(plot1,'LineWidth',2)
xlabel('Freqency range (Hz)'); ylabel('Transmittance (dB)');
title('Swept Sine results for displacement')
grid on
legend 'mass1' 'mass2' 'last mass1' 'last mass2'
set(gca,'fontsize',20) 
%% Non-dimensionalised Results
% non-dimensionalise the frequency 
NDf=swept_sine_range/w2;
% Plot the time series
figure
plot1=plot(NDf,mag2db(m1_disp),NDf,mag2db(m2_disp)...
    ,NDf,mag2db(lm1_disp),NDf,mag2db(lm2_disp));
set(plot1,'LineWidth',2)
xlabel('\eta'); ylabel('Transmittance (dB)');
title('Transmittance curve for weakly NL MIM model with NL = 1600')
grid on
legend 'mass1' 'mass2' 'last mass1' 'last mass2'
set(gca,'fontsize',24) 

%% Transmission
% U=Xn/X1
U=abs(lm2_disp)./(m2_disp);
figure
plot(swept_sine_range,20*log10((U)),'b')
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
function dxdt=rhs(t,x,omega,k3)
        mass1=0.1;		% [kg]
        mass2=0.5*mass1;
        stiff1=1000;    % [N/m]
        stiff2=1500;
        stiff3=k3;
        damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.002;
        f=1; %*(stepfun(t,0)-stepfun(t,0.01));
        w=omega; % Hz, forcing frequency 
        %----1 unit cell-----
        u1=x(1);    %disp mass1
        du1=x(2);    %velo mass1
        v1=x(3);   %disp mass2
        dv1=x(4);  % velo mass2
        %----2 unit cell-----
        u2=x(5);    %disp mass1
        du2=x(6);    %velo mass1
        v2=x(7);   %disp mass2
        dv2=x(8);  % velo mass2
        %----3 unit cell-----
        u3=x(9);    %disp mass1
        du3=x(10);    %velo mass1
        v3=x(11);   %disp mass2
        dv3=x(12);  % velo mass2
        %----4 unit cell-----
        u4=x(13);    %disp mass1
        du4=x(14);    %velo mass1
        v4=x(15);   %disp mass2
        dv4=x(16);  % velo mass2
        %----5 unit cell-----
        u5=x(17);    %disp mass1
        du5=x(18);    %velo mass1
        v5=x(19);   %disp mass2
        dv5=x(20);  % velo mass2
     
        %---------------------------------------
        %---------------1 unit cell
        % first mass
        dxdt_1 = du1;
        dxdt_2 = -((2*damp1+damp2)/mass1)*du1- ((2*stiff1)/mass1)*u1-(stiff2/mass1)*(u1-v1) -...
            (stiff3/mass1)*(u1-v1)^3+(damp2/mass1)*dv1 + (stiff1/mass1)*u2 + (damp1/mass1)*du2 ...
            +(f/mass1)*sin(2*pi*w*t);
        % second mass
        dydt_1= dv1;
        dydt_2= -(stiff2/mass2)*(v1-u1)-(stiff3/mass2)*(v1-u1)^3 - (damp2/mass2)*dv1 + (damp2/mass2)*du1;
        %---------------------------------------
        %---------------2 unit cell 
        % first mass
        dxdt_3 = du2;
        dxdt_4 = -((2*damp1+damp2)/mass1)*du2- ((2*stiff1)/mass1)*u2-(stiff2/mass1)*(u2-v2) -...
            (stiff3/mass1)*(u2-v2)^3+(damp2/mass1)*dv2 + (stiff1/mass1)*u1 + (damp1/mass1)*du1...
            + (stiff1/mass1)*u3 + (damp1/mass1)*du3;
        % second mass
        dydt_3= dv2;
        dydt_4= -(stiff2/mass2)*(v2-u2)-(stiff3/mass2)*(v2-u2)^3 - (damp2/mass2)*dv2 + (damp2/mass2)*du2;
        %---------------------------------------
        %---------------3 unit cell 
        % first mass
        dxdt_5 = du3;
        dxdt_6 = -((2*damp1+damp2)/mass1)*du3- ((2*stiff1)/mass1)*u3-(stiff2/mass1)*(u3-v3) -...
            (stiff3/mass1)*(u3-v3)^3+(damp2/mass1)*dv3 + (stiff1/mass1)*u2 + (damp1/mass1)*du2...
            + (stiff1/mass1)*u4 + (damp1/mass1)*du4;
        % second mass
        dydt_5= dv3;
        dydt_6= -(stiff2/mass2)*(v3-u3)-(stiff3/mass2)*(v3-u3)^3 - (damp2/mass2)*dv3 + (damp2/mass2)*du3;
        %---------------------------------------
        %---------------4 unit cell 
        % first mass
        dxdt_7 = du4;
        dxdt_8 = -((2*damp1+damp2)/mass1)*du4- ((2*stiff1)/mass1)*u4-(stiff2/mass1)*(u4-v4) -...
            (stiff3/mass1)*(u4-v4)^3+(damp2/mass1)*dv4 + (stiff1/mass1)*u3 + (damp1/mass1)*du3...
            + (stiff1/mass1)*u5 + (damp1/mass1)*du5;
        % second mass
        dydt_7= dv4;
        dydt_8= -(stiff2/mass2)*(v4-u4)-(stiff3/mass2)*(v4-u4)^3 - (damp2/mass2)*dv4 + (damp2/mass2)*du4;
        %---------------------------------------
        %---------------5 unit cell (last cell)
        % first mass
        dxdt_9 = du5;
        dxdt_10 = -((2*damp1+damp2)/mass1)*du5- ((2*stiff1)/mass1)*u5-(stiff2/mass1)*(u5-v5) -...
            (stiff3/mass1)*(u5-v5)^3+(damp2/mass1)*dv5 + (stiff1/mass1)*u4 + (damp1/mass1)*du4;
        % second mass
        dydt_9= dv5;
        dydt_10= -(stiff2/mass2)*(v5-u5)-(stiff3/mass2)*(v5-u5)^3 - (damp2/mass2)*dv5 + (damp2/mass2)*du5;
        %---------------------------------------
        % -----------Final Solution-------------
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2;dxdt_3;dxdt_4;dydt_3;dydt_4;...
            dxdt_5; dxdt_6; dydt_5; dydt_6;dxdt_7;dxdt_8;dydt_7;dydt_8;...
            dxdt_9; dxdt_10; dydt_9; dydt_10];
end