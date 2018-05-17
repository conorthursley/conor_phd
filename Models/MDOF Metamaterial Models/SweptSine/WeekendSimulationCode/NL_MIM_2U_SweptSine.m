% Plots the transient response of a forced 2 unit cell MDOF metamaterial model

tic
% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;    % [s] delta t
% for loop parameters
t_end=1000;   % t limit
t=0:dt:t_end;      % [s] time scaled/dx of sined
t_find=600; % the time to safely assume SS has been reached 600 seconds after initial transient begins
p=find(t==600); q=find(t==t_end);

mass1=0.1;		% [kg]
mass2=0.5*mass1;
stiff1=1000;    % [N/m]
stiff2=1500;

w2=sqrt(stiff2/mass2)/(2*pi);
theta=mass2/mass1;


% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
z=zeros(1,2*4); % n=2 and there are 4 DOF per unit cell, hence 4 initial conditions;
% Set the frequency range
freq_step=0.005;
swept_sine_range=11.58:freq_step:14.34; % range from 10 Hz to 45 Hz in steps of 0.25 Hz

% set the nonlinear strength
sigma=[0 100 200 400 800 1600]*stiff2;
% set the storage cell
amplitude=zeros(length(swept_sine_range),length(z)); %cell to store amps of displacement and velocity

% Solve the model
for j=1:length(sigma)
    k3=sigma(j);
    amplitude(j,:,:)=k3;
    parfor i=1:length(swept_sine_range)
        t=0:dt:t_end;
        omega=swept_sine_range(i);
        options=odeset('InitialStep',dt,'MaxStep',dt);
        [t,result]=ode45(@(t,z) rhs(t,z,omega,k3),t,z,options);
        
        %change result to show the steady state portion of the time history
        x=result(p:q,:); % x becomes the steady state result
        amplitude(i,:)=max(x);
    end
    t_new=t(p:q);
    toc
    
    % Results
    % extract displacement amplitudes from vector, "amplitude"
    % amplitude = [displacement1 velo1 disp2 velo2 disp'n' velo'n' .....]
    % ----------First unit cell--------------
    %  displacement
    m1_disp=amplitude(:,1);
    m2_disp=amplitude(:,3);
    %  velocity
    m1_velo=amplitude(:,2);
    m2_velo=amplitude(:,4);
    %   ----------Last unit cell---------------
    %  displacement, 'last mass' - 'lm'
    lm1_disp=amplitude(:,end-3);
    lm2_disp=amplitude(:,end-1);
%     velocity
    lm1_velo=amplitude(:,end-2);
    lm2_velo=amplitude(:,end);
    
    % Non-dimensionalised Results
    % non-dimensionalise the frequency
    NDf=swept_sine_range/w2;
    %Plot the time series
    if j==1
        figure
    else
        hold on
    end
    
    plot1=plot(NDf,mag2db(lm1_disp),NDf,mag2db(lm2_disp));
    set(plot1,'LineWidth',2)
    xlabel('\eta'); ylabel('Transmittance (dB)');
    title('Transmittance curve for 2U system MIM model 10 amp')
    grid on
    legend 'mass1' 'mass2' 'last mass1' 'last mass2'
    set(gca,'fontsize',20)
    hold on
end
% Save the results
save 2UweaklyNLamp10_modeGroup1.mat amplitude swept_sine_range t t_new
savefig('2UweaklyNLamp10_modeGroup1.fig')
% Mass-Spring-Damper system
%The equations for the mass spring damper system have to be defined
%separately so that the ODE45 solver can call it.
function dxdt=rhs(t,x,omega,k3)
mass1=0.1;		% [kg]
mass2=0.5*mass1;
stiff1=1000;    % [N/m]
stiff2=1500;
stiff3=k3;
damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
damp2=0.002;
f=10; %*(stepfun(t,0)-stepfun(t,0.01));
w=omega; % Hz, forcing frequency
%----first unit cell-----
u1=x(1);    %disp mass1
du1=x(2);    %velo mass1
v1=x(3);   %disp mass2
dv1=x(4);  % velo mass2
%----second unit cell-----
u2=x(5);    %disp mass1
du2=x(6);    %velo mass1
v2=x(7);   %disp mass2
dv2=x(8);  % velo mass2

%---------------------------------------
%first unit cell
%first mass
dxdt_1 = du1;
dxdt_2 = -((2*damp1+damp2)/mass1)*du1- ((2*stiff1)/mass1)*u1-(stiff2/mass1)*(u1-v1) -...
    (stiff3/mass1)*(u1-v1)^3+(damp2/mass1)*dv1 + (stiff1/mass1)*u2 + (damp1/mass1)*du2 ...
    +(f/mass1)*sin(2*pi*w*t);
%second mass
dydt_1= dv1;
dydt_2= -(stiff2/mass2)*(v1-u1)-(stiff3/mass2)*(v1-u1)^3 - (damp2/mass2)*dv1 + (damp2/mass2)*du1;
%---------------------------------------
%second unit cell (last cell)
%first mass
dxdt_3 = du2;
dxdt_4 = -((2*damp1+damp2)/mass1)*du2- ((2*stiff1)/mass1)*u2-(stiff2/mass1)*(u2-v2) -...
    (stiff3/mass1)*(u2-v2)^3+(damp2/mass1)*dv2 + (stiff1/mass1)*u1 + (damp1/mass1)*du1;
%second mass
dydt_3= dv2;
dydt_4= -(stiff2/mass2)*(v2-u2)-(stiff3/mass2)*(v2-u2)^3 - (damp2/mass2)*dv2 + (damp2/mass2)*du2;
%-----------Final Solution-------------
dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2;dxdt_3;dxdt_4;dydt_3;dydt_4];
end