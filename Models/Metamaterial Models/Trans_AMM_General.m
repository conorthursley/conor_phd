%% Plots the transient response of a forced 2DOF mass spring damper system 

%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:10;      % [s] time scale
%% Matrices
M=[0.1 0;0 0.05];
K=[250 -100;-100 100];
omega=logspace(-1,2,length(t));

f=[1; 1];

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0;
initial_dxdt = 0;
initial_y    = 0;
initial_dydt = 0;
initial_x1    = 0;
initial_dx1dt = 0;
initial_y1    = 0;
initial_dy1dt = 0;
z=[initial_x initial_dxdt initial_y initial_dydt initial_x1 initial_dx1dt initial_y1 initial_dy1dt];
%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45(@rhs, t, z, options);
%% Import comparison 
% M='U:\_PhD\Datathief\MDOF_freeResponse-InmanEngVib\figure12_mode1.csv';
% data=csvread(M,1,0);
% N='U:\_PhD\Datathief\MDOF_freeResponse-InmanEngVib\figure12_mode2.csv';
% data1=csvread(N,1,0);
%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),t,x(:,3));
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' 'Datathief mode1' 'Datathief mode2'
set(gca,'fontsize',20) 
%%% Calculate the PSD of the time series
FFTsize=1024;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x,hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');
%% Amplitude and resonance 
Xamp=zeros(length(omega),2)';
Wn=zeros(length(omega),2)';
for jj=1:length(omega)
    Wn(:,jj)=[omega(jj);omega(jj)];
    [amp]=forced_vibration(K,M,f,Wn(jj));
    Xamp(:,jj)=amp;
end
figure
plot(omega,abs(Xamp))
%% TF estimate
% takes in input and output signal
% bandpass filter input signal
wind = kaiser(length(x(:,7)),105);
signal=x(:,1);
% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(signal,x(:,7),[],[],[],1/(dt));
% plot
figure
% hold on
graph=plot(frequencies,(20*log10(abs(txy))),'r');
% axis([0 12 -Inf Inf])
grid on
% title('Transfer function of metamaterial configurations','FontSize',14)
xlabel('Frequency, \omega','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',14)
% axis([0 10 -Inf 0])
legend({'linear AMM','noise insulation panel','nonlinear case 1'},'FontSize',14)
set(graph,'LineWidth',1.5);
%% semilog txy
figure
% hold on
graph1=semilogx(frequencies,abs(txy),'b');
% axis([0 3 0 Inf])
grid on
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',20)
legend({'linear AMM numerical result','nonlinear case 1 AMM numerical result','phononic crystal numerical result','nonlinear case 2 AMM numerical result'})
set(graph1,'LineWidth',2);
%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=rhs(t,x)
        mass1=0.1;		% [kg]
        mass2=0.3*mass1;
        stiff1=1e3;    % [N/m]
        stiff2=0.1*stiff1;
        damp=0;     % [Ns/m] keep as a small number to fix solver errors
        f=1; %*(stepfun(t,0)-stepfun(t,0.01));
     
        u1=x(1);
        u1_v=x(2);
        v1=x(3);
        v1_v=x(4);
        u2=x(5);
        u2_v=x(6);
        v2=x(7);
        v2_v=x(8);
        
        dxdt_1 = u1_v;
        dxdt_2 = - ((2*stiff1+stiff2)/mass1)*u1 +(stiff2/mass1)*v1 + (stiff1/mass1)*u2+(f/mass1)*cos(2*pi*10*t);
        dydt_1= v1_v;
        dydt_2= (stiff2/mass2)*u1 - (stiff2/mass2)*v1;
        dxdt_3 = u2_v;
        dxdt_4 = - ((stiff1+stiff2)/mass1)*u2 +(stiff2/mass1)*v2 + (stiff1/mass1)*u1;
        dydt_3= v2_v;
        dydt_4= (stiff2/mass2)*u2 - (stiff2/mass2)*v2;

        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2;dxdt_3;dxdt_4;dydt_3;dydt_4];
end
%% forced vibration 
function X = forced_vibration(K,M,f,omega)
% Function to calculate steady state amplitude of
% a forced linear system.
% K is nxn the stiffness matrix
% M is the nxn mass matrix
% f is the n dimensional force vector
% omega is the forcing frequency, in radians/sec.
% The function computes a vector X, giving the amplitude of
% each degree of freedom
%
X = (K-M*omega^2)\f;

end