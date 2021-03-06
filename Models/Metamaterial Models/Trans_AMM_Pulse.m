%% Plots the transient response of a forced 2DOF mass spring damper system 

%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:30;      % [s] time scale
omega = logspace(-1,3, length(t));
%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0;
initial_dxdt = 0;
initial_y    = 0;
initial_dydt = 0;
z=[initial_x initial_y initial_dxdt initial_dydt];
%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt,'RelTol',1e-4);
[t,x]=ode45(@rhs, t, z, options);

%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),'r*',t,x(:,2),'b+');
set(plot1,'LineWidth',1.5)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' %'Datathief mode1' 'Datathief mode2'
set(gca,'fontsize',20) 
%%% Calculate the PSD of the time series
FFTsize=1024;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x,hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');

%% Amplitude Response
% parameters should be passed through ODE45 to avoid reptitive changing of
% inputs
mass1=101.1e-3;		% [kg]
mass2=45.47e-3;
stiff1=117;    % [N/m]
stiff2=74;
damp1=0.0002;     % [Ns/m] keep as a small number to fix solver errors
damp2=0.0002;
w=27*2*pi; % driving frequency
% matrix allocation
M=[mass1 0; 0 mass2];   % mass matrix
C=[damp1+damp2 -damp2;-damp2 damp2];  % Damping matrix
K=[stiff1+stiff2 -stiff2;-stiff2 stiff2];
% matrix in the form seen here
% http://www.brown.edu/Departments/Engineering/Courses/En4/Notes/vibrations_mdof/vibrations_mdof.htm
N=[eye(2) zeros(2);zeros(2) M];
D=[zeros(2) -eye(2);K C];
B=[1;0];
f=[0;0;B];
for jj=1:length(omega)
[amp,phase] = damped_forced_vibration(D,N,f,omega(jj));
X(:,jj)=amp;
end
figure
loglog(omega,X)
%% FFT
figure
dt=abs(mean(diff(t)));  %average time step done 
Fs=1/(dt);
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(t); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(x(:,1),m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft1); 
f=Fs*(0:(n/2))/n;
freq1=fr(1:floor(m/2));
P1=fourier(1:floor(m/2));
plot(freq1,P1)
% plot(flop,abs(y),'LineWidth',2)
title('FFT')
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',20)

%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x)
        mass1=0.1;		% [kg]
        mass2=0.5*mass1;
        stiff1=1000;    % [N/m]
        stiff2=1500;
        damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.002;
        w=10*2*pi; % driving frequency
        % matrix allocation
        M=[mass1 0; 0 mass2];   % mass matrix
        C=[damp1+damp2 -damp2;-damp2 damp2];  % Damping matrix
        K=[stiff1+stiff2 -stiff2;-stiff2 stiff2]; %stiffness matrix
        B=[1;0]; % forcing input vector
        t1=1; t2=1.1; % pulse duration
        
        
        
        % state space formulation
        A1=[zeros(2) eye(2); -inv(M)*K -inv(M)*C];
        f=inv(M)*B;
        
        % solution 
        dxdt=A1*x+[0;0;f]*sin(w*t); %(stepfun(t,t1)-stepfun(t,t2));
    end
%% Forced Vibraiton Amplitude function 
function [amp,phase] = damped_forced_vibration(D,M,f,omega)

% Function to calculate steady state amplitude of
% a forced linear system.
% D is 2nx2n the stiffness/damping matrix
% M is the 2nx2n mass matrix
% f is the 2n dimensional force vector
% omega is the forcing frequency, in radians/sec.
% The function computes a vector �amp�, giving the amplitude of
% each degree of freedom, and a second vector �phase�,
% which gives the phase of each degree of freedom

Y0 = (D+M.*1i*omega)\f;  % The i here is sqrt(-1)
% We dont need to calculate Y0bar - we can just change the sign of
% the imaginary part of Y0 using the 'conj' command
for j =1:length(f)/2
    amp(j) = sqrt(Y0(j)*conj(Y0(j)));
    phase(j) = log(conj(Y0(j))/Y0(j))/(2*i);
end

end