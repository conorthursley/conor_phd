%% Plots the transient response of a forced mass spring damper system 
clear all
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:100;      % [s] time scale
tic
%% Initial conditions: x(0) = 0, x'(0)=0 
initial_x    = 0;
initial_dxdt = 0;
w=logspace(-1, 0.5, 200);
w_rev=fliplr(w);
vec=zeros(length(w),1)';
vec_rev=vec;
%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
parfor i =1:length(w)
    w_f=w(i);
    t=0:dt:100;
    [t,x]=ode45(@(t,y) rhs(t,y,w_f), t, [initial_x initial_dxdt],options );
    rms=sqrt(mean(x(:,1).^2));
    vec(i)=rms;
    w_f=w_rev(i);
    [t,x]=ode45(@(t,y) rhs(t,y,w_f), t, [initial_x initial_dxdt],options );
    rms=sqrt(mean(x(:,1).^2));
    vec_rev(i)=rms;
end

toc
%% Plot the results
figure
plot(w,vec,w,vec_rev)
grid on

%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x,w)
        mass1=1;		% [kg]
        stiff1=1;    % [N/m]
        stiff2=0;
        damp=0.32;     % [Ns/m] keep as a small number to fix solver errors
        f=2.5; %*sin(2*pi*10*t);            % [N] amplitude of driving force
        
        
        dxdt_1 = x(2);
        dxdt_2 = -(damp/(mass1))*x(2) - (stiff1/mass1)*x(1) -(stiff2/mass1).*x(1)^3 + (f/mass1)*cos(w*t);

        dxdt=[dxdt_1; dxdt_2];
    end