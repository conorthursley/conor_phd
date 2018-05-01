%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t_end=4000;   % t limit
t=0:dt:t_end;      % [s] time scale

omega=25;
%work period/cycle
%from 8000 seconds
t_b=3500;
t_cycle=2*pi*sqrt(1/((omega*(2*pi)))^2);
t_a=((t_b+t_cycle));
[m_min,i_min]=min(abs(t(:)-t_a));
p=i_min; q=find(t==t_b); 

%% Initial conditions: x(0) = 0, x'(0)=0
initial_x    = 0;
initial_dxdt = 0;
K=1000;
mass1=0.1;
force=1;
damp=0.002;
z=[initial_x initial_dxdt];
w_nat=sqrt(K/mass1);

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45(@(t,z) rhs(t,z,omega),t,z,options);
t_new=t(q:p);
toc

%% Work and Energy functions
%-----------Results----------------------
% extract steady state values from vector, "x"
% x = [displacement1 velo1]
x_new=x(q:p,:);
% displacement
m1_disp=x_new(:,1);
% velocity 
m1_velo=x_new(:,2);

%-----------Sum of Forces-------------------
% harmonic solution assumed for F and x_disp
% Input Force
In=(force)*sin(omega*2*pi*t_new);
% Output motion 
Out=-mass1*((omega*2*pi)^2)*(m1_disp) + K*(m1_disp);
%plot
figure
plot(t_new,In,'b',t_new,Out,'r--')

%------------Work-----------------------------
% Looking at the dynamic model, apply the energy approach for KE and PE and
% find the total KE and PE of the system
% -------------KE------------
% Kinetic Energy = 0.5*m_i*v_i^2
KE=0.5*mass1*((m1_velo).^2);
% -------------PE------------
% Potential Energy = 0.5*k_i*u_i^2
PE=0.5*K*((m1_disp).^2);

% total mechanical energy
ME=-PE+KE;
workIn=(force*sin(omega*2*pi*t_new)).*((m1_disp));
%plot
figure
plot(t_new,workIn,'b',t_new,ME,'r--')


%% Work and Energy Figure
figure
plot1=plot(t_new,KE,'r',t_new,PE,'b',t_new,ME,'k',[t_b,t_a],[y1, y1],'m:',...
    t_new,workIn,'g',[t_b,t_a],[y2, y2],'c:');
xlabel('t'); ylabel('Work (J)');
set(plot1,'LineWidth',2.5)
title(['Work/Energy Calculations at ', num2str(omega), ' Hz'])
grid on
axis([t_b t_a -Inf Inf]);
legend 'KE' 'PE' 'Work Output' 'Work Output RMS' 'Work Input'  'Work Input RMS'
set(gca,'fontsize',24) 
%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x,omega)
        mass1=0.1;		% [kg]
        stiff1=1000;    % [N/m]
        damp=0.002;     % [Ns/m] keep as a small number to fix solver errors
        f=1;           % [N] amplitude of driving force
        w=omega;

        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - (stiff1/mass1)*x(1) + (f/mass1)*sin((w*2*pi)*t);

        dxdt=[dxdt_1; dxdt_2];
    end