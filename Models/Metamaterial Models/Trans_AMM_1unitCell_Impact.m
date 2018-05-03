%% 1U Impact metamaterial model
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
z=zeros(1,1*4); % n=1 and there are 4 DOF per unit cell, hence 4 initial conditions;

%% set the frequency
omega=20;   % Frequency [Hz]
%% Solve the model

options=odeset('InitialStep',dt,'MaxStep',dt);
[t,result]=ode45(@(t,z) rhs(t,z,omega),t,z,options);

% change result to show the steady state portion of the time history
x=result(p:q,:); % x becomes the steady state result

t_new=t(p:q);
toc

%% Plot the results
% Plot the time series
figure
plot1=plot(t_new,x);
set(plot1,'LineWidth',2)
xlabel('Time,s'); ylabel('Displacement, Velocity');
title('Impact Oscillator MIM model')
grid on
set(gca,'fontsize',20) 
% Plot the phase plane
figure
plot1=plot(x(:,1),x(:,2));
set(plot1,'LineWidth',2)
xlabel('Displacement'); ylabel('Velocity');
title('Impact Oscillator MIM model Phase Portrait')
grid on
set(gca,'fontsize',20) 
%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=rhs(t,x,omega)
        mass1=0.1;		% [kg]
        mass2=0.5*mass1;
        stiff1=1000;    % [N/m]
        stiff2=1500;
        damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.002;
        f=1; %*(stepfun(t,0)-stepfun(t,0.01));
        w=omega; % Hz, forcing frequency 
        
        
        %----first unit cell-----
        u1=x(1);    %disp mass1
        du1=x(2);    %velo mass1
        v1=x(3);   %disp mass2
        dv1=x(4);  % velo mass2
        %---------------------------------------
        % first unit cell
        % first mass
        dxdt_1 = du1;
        dxdt_2 = -((2*damp1+damp2)/mass1)*du1- ((2*stiff1+stiff2)/mass1)*u1...
           +(damp2/mass1)*dv1 + (stiff2/mass1)*v1 + (f/mass1)*sin(2*pi*w*t);
        % second mass
        dydt_1= dv1;
        dydt_2= -(stiff2/mass2)*(v1) - (damp2/mass2)*dv1 + (damp2/mass2)*du1...
            + (stiff2/mass2)*(u1);
        %---------------------------------------
        % -----------Final Solution-------------
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2];
end