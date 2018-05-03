%% [5 CELL]   Plots the transient response of a forced 5 unit cell MDOF metamaterial model
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
z=zeros(1,10*4); % n=5 and there are 4 DOF per unit cell, hence 4 initial conditions;
%% Set the frequency range
freq_step=0.25;
swept_sine_range=freq_step:freq_step:60; % range from 10 Hz to 45 Hz in steps of 0.25 Hz

%% set the nonlinear strength
sigma=[0 100 200 400 800 1600]*stiff2;
%% set the storage cell
amplitude=zeros(length(swept_sine_range),length(z)); %cell to store amps of displacement and velocity

%% Solve the model
for j=1:length(sigma)
    k3=sigma(j);
%     amplitude(j,:,:)=k3;
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
    
%% Non-dimensionalised Results
    % non-dimensionalise the frequency
    NDf=swept_sine_range/w2;
    % Plot the time series
    %create figure
    if j==1
        figure
    else
        hold on
    end
    % initiate plot
    plot1=plot(NDf,mag2db(lm1_disp),NDf,mag2db(lm2_disp));
    set(plot1,'LineWidth',2)
    xlabel('\eta'); ylabel('Transmittance (dB)');
    title('Transmittance curve for 10U system MIM model')
    grid on
    legend 'mass1' 'mass2' 'last mass1' 'last mass2'
    set(gca,'fontsize',20)
    hold on
end
%% Save the results
save 10UweaklyNL.mat amplitude swept_sine_range t t_new
savefig('10UweaklyNL.fig')
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
        %----6 unit cell-----
        u6=x(21);    %disp mass1
        du6=x(22);    %velo mass1
        v6=x(23);   %disp mass2
        dv6=x(24);  % velo mass2
        %----7 unit cell-----
        u7=x(25);    %disp mass1
        du7=x(26);    %velo mass1
        v7=x(27);   %disp mass2
        dv7=x(28);  % velo mass2
        %----8 unit cell-----
        u8=x(29);    %disp mass1
        du8=x(30);    %velo mass1
        v8=x(31);   %disp mass2
        dv8=x(32);  % velo mass2
        %----9 unit cell-----
        u9=x(33);    %disp mass1
        du9=x(34);    %velo mass1
        v9=x(35);   %disp mass2
        dv9=x(36);  % velo mass2
        %----10 unit cell-----
        u10=x(37);    %disp mass1
        du10=x(38);    %velo mass1
        v10=x(39);   %disp mass2
        dv10=x(40);  % velo mass2
     
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
         %---------------5 unit cell 
        % first mass
        dxdt_9 = du5;
        dxdt_10 = -((2*damp1+damp2)/mass1)*du5- ((2*stiff1)/mass1)*u5-(stiff2/mass1)*(u5-v5) -...
            (stiff3/mass1)*(u5-v5)^3+(damp2/mass1)*dv5 + (stiff1/mass1)*u4 + (damp1/mass1)*du4...
            + (stiff1/mass1)*u6 + (damp1/mass1)*du6;
        % second mass
        dydt_9= dv5;
        dydt_10= -(stiff2/mass2)*(v5-u5)-(stiff3/mass2)*(v5-u5)^3 - (damp2/mass2)*dv5 + (damp2/mass2)*du5;
        %---------------------------------------
        %---------------6 unit cell 
        % first mass
        dxdt_11 = du6;
        dxdt_12 = -((2*damp1+damp2)/mass1)*du6- ((2*stiff1)/mass1)*u6-(stiff2/mass1)*(u6-v6) -...
            (stiff3/mass1)*(u6-v6)^3+(damp2/mass1)*dv6 + (stiff1/mass1)*u5 + (damp1/mass1)*du5...
            + (stiff1/mass1)*u7 + (damp1/mass1)*du7;
        % second mass
        dydt_11= dv6;
        dydt_12= -(stiff2/mass2)*(v6-u6)-(stiff3/mass2)*(v6-u6)^3 - (damp2/mass2)*dv6 + (damp2/mass2)*du6;
        %---------------------------------------
        %---------------7 unit cell 
        % first mass
        dxdt_13 = du7;
        dxdt_14 = -((2*damp1+damp2)/mass1)*du7- ((2*stiff1)/mass1)*u7-(stiff2/mass1)*(u7-v7) -...
            (stiff3/mass1)*(u7-v7)^3+(damp2/mass1)*dv7 + (stiff1/mass1)*u6 + (damp1/mass1)*du6...
            + (stiff1/mass1)*u8 + (damp1/mass1)*du8;
        % second mass
        dydt_13= dv7;
        dydt_14= -(stiff2/mass2)*(v7-u7)-(stiff3/mass2)*(v7-u7)^3 - (damp2/mass2)*dv7 + (damp2/mass2)*du7;
        %---------------------------------------
        %---------------8 unit cell 
        % first mass
        dxdt_15 = du8;
        dxdt_16 = -((2*damp1+damp2)/mass1)*du8- ((2*stiff1)/mass1)*u8-(stiff2/mass1)*(u8-v8) -...
            (stiff3/mass1)*(u8-v8)^3+(damp2/mass1)*dv8 + (stiff1/mass1)*u7 + (damp1/mass1)*du7...
            + (stiff1/mass1)*u9 + (damp1/mass1)*du9;
        % second mass
        dydt_15= dv8;
        dydt_16= -(stiff2/mass2)*(v8-u8)-(stiff3/mass2)*(v8-u8)^3 - (damp2/mass2)*dv8 + (damp2/mass2)*du8;
        %---------------------------------------
         %---------------9 unit cell 
        % first mass
        dxdt_17 = du9;
        dxdt_18 = -((2*damp1+damp2)/mass1)*du9- ((2*stiff1)/mass1)*u9-(stiff2/mass1)*(u9-v9) -...
            (stiff3/mass1)*(u9-v9)^3+(damp2/mass1)*dv9 + (stiff1/mass1)*u8 + (damp1/mass1)*du8...
            + (stiff1/mass1)*u10 + (damp1/mass1)*du10;
        % second mass
        dydt_17= dv9;
        dydt_18= -(stiff2/mass2)*(v9-u9)-(stiff3/mass2)*(v9-u9)^3 - (damp2/mass2)*dv9 + (damp2/mass2)*du9;
        %---------------------------------------
        %---------------10 unit cell (last cell)
        % first mass
        dxdt_19 = du10;
        dxdt_20 = -((2*damp1+damp2)/mass1)*du10- ((2*stiff1)/mass1)*u10-(stiff2/mass1)*(u10-v10) -...
            (stiff3/mass1)*(u10-v10)^3+(damp2/mass1)*dv10 + (stiff1/mass1)*u9 + (damp1/mass1)*du9;
        % second mass
        dydt_19= dv10;
        dydt_20= -(stiff2/mass2)*(v10-u10)-(stiff3/mass2)*(v10-u10)^3 - (damp2/mass2)*dv10 + (damp2/mass2)*du10;
        %---------------------------------------
        % -----------Final Solution-------------
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2;dxdt_3;dxdt_4;dydt_3;dydt_4;...
            dxdt_5; dxdt_6; dydt_5; dydt_6;dxdt_7;dxdt_8;dydt_7;dydt_8;...
            dxdt_9; dxdt_10; dydt_9; dydt_10;dxdt_11; dxdt_12; dydt_11; ...
            dydt_12;dxdt_13;dxdt_14;dydt_13;dydt_14;...
            dxdt_15; dxdt_16; dydt_15; dydt_16;dxdt_17;dxdt_18;dydt_17;dydt_18;...
            dxdt_19; dxdt_20; dydt_19; dydt_20];
end