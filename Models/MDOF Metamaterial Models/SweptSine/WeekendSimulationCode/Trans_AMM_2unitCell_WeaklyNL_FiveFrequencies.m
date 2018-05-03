%% Plots the transient response of a forced 2DOF mass spring damper system 

tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;    % [s] delta t
% for loop parameters
t_end=1500;   % t limit
t=0:dt:t_end;      % [s] time scale
t_find=1200; % the time to safely assume SS has been reached 600 seconds after initial transient begins
p=find(t==t_find); q=find(t==t_end);

mass1=0.1;		% [kg]
mass2=mass1*0.5;
stiff1=1000;    % [N/m]
stiff2=1.5*stiff1;
force=1;
w2=sqrt(stiff2/mass2)/(2*pi);
w0=sqrt(stiff2/mass2);
theta=mass2/mass1;

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
z=zeros(1,2*4); % n=2 and there are 4 DOF per unit cell, hence 4 initial conditions;
%% Set the frequency range
f1=10;
f2=17;
f3=28;
f4=37;
f5=43;

freq_range=[f1 f2 f3 f4 f5]; % range from 10 Hz to 43 Hz 
% select which mass we want to observe the frequency response for 
res=1; % result of the m1 displacement = 1,
% m2 displacement = 3
% freq_results=zeros(length(p:q),length(freq_range));% p and q are the new time start/end to cut 
%% set the nonlinear strength vector 
stiff3=[0 100 200 400 800 1600]*stiff2;
%% Solve the model
for jjj=1:length(stiff3)
    k3=stiff3(jjj);
    parfor i=1:length(freq_range)
        omega=freq_range(i);
        t=0:dt:t_end;      % [s] time scale
        options=odeset('InitialStep',dt,'MaxStep',dt);
        [t,result]=ode45(@(t,z) rhs(t,z,omega,k3),t,z,options);
        
        % change result to show the steady state portion of the time history
        x=result;%(p:q,:); % x becomes the steady state result
        freq_results1(:,i)=x(:,1); %store the displacement history of 1 cell mass1
        freq_results2(:,i)=x(:,5); %store the displacement history of 2 cell mass1
        freq_results3(:,i)=x(:,3); %store the displacement history of 1 cell mass2
    end
    
    toc
    %% Results
    % results are stored in freq_results vector
    % need to normalise with respect to the static displacement
    % static displacement = disp/(F/k1)=disp/(1/1000)
    for ii=1:length(freq_range)
        disp1(:,ii)=freq_results1(:,ii)/(force/stiff1);
        disp2(:,ii)=freq_results2(:,ii)/(force/stiff1);
        disp3(:,ii)=freq_results3(:,ii)/(force/stiff1);
    end
    % subplots of the five frequencies together
    R=p+ceil(25*w0);
    % first two masses
    figure
    for ii=1:length(freq_range)
        subplot(length(freq_range),1,ii)
        plot(t(p:R),disp1(p:R,ii),'b',t(p:R),disp2(p:R,ii),'r--');
        str=(['\eta = ',num2str(freq_range(ii)/w2)]);
        title(str)
    end
    % second mass of the first cell
    figure
    for ii=1:length(freq_range)
        subplot(length(freq_range),1,ii)
        plot(t(p:R),disp3(p:R,ii),'m');
        str=(['\eta = ',num2str(freq_range(ii)/w2)]);
        title(str)
    end
    
    % both together
    h=figure;
    %counters for subplots
    j=1;
    jj=2;
    for ii=1:length(freq_range)
        % first masses
        subplot(length(freq_range),2,j)
        plot(t(p:R),disp1(p:R,ii),'b',t(p:R),disp2(p:R,ii),'r--');
        str=(['\eta = ',num2str(freq_range(ii)/w2)]);
        title(str)
        % second mass
        subplot(length(freq_range),2,jj)
        plot(t(p:R),disp3(p:R,ii),'m');
        str=(['\eta = ',num2str(freq_range(ii)/w2)]);
        title(str)
        j=j+2;
        jj=jj+2;
    end
    fileFig=['2U_5Freq_',num2str(k3/stiff2),'.fig'];
    savefig(h,fileFig)
    filestr=['2U_5Freq_',num2str(k3/stiff2),'.mat'];
    save(filestr, 'disp1', 'disp2', 'disp3', 'R')
end

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
        % first unit cell
        % first mass
        dxdt_1 = du1;
        dxdt_2 = -((2*damp1+damp2)/mass1)*du1- ((2*stiff1)/mass1)*u1-(stiff2/mass1)*(u1-v1) -...
            (stiff3/mass1)*(u1-v1)^3+(damp2/mass1)*dv1 + (stiff1/mass1)*u2 + (damp1/mass1)*du2 ...
            +(f/mass1)*sin(2*pi*w*t);
        % second mass
        dydt_1= dv1;
        dydt_2= -(stiff2/mass2)*(v1-u1)-(stiff3/mass2)*(v1-u1)^3 - (damp2/mass2)*dv1 + (damp2/mass2)*du1;
        %---------------------------------------
        % second unit cell (last cell)
        % first mass
        dxdt_3 = du2;
        dxdt_4 = -((2*damp1+damp2)/mass1)*du2- ((2*stiff1)/mass1)*u2-(stiff2/mass1)*(u2-v2) -...
            (stiff3/mass1)*(u2-v2)^3+(damp2/mass1)*dv2 + (stiff1/mass1)*u1 + (damp1/mass1)*du1;
        % second mass
        dydt_3= dv2;
        dydt_4= -(stiff2/mass2)*(v2-u2)-(stiff3/mass2)*(v2-u2)^3 - (damp2/mass2)*dv2 + (damp2/mass2)*du2;
        % -----------Final Solution-------------
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2;dxdt_3;dxdt_4;dydt_3;dydt_4];
end